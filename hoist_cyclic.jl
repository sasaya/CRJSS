using JuMP, HiGHS, GLMakie #CairoMakie
using JSON

function main(W,R,N,Capacity,V,L,U,E,D;num_of_hoist = 1:2)
    function Qlamda(y) 
        return getindex.(findall(x -> x == y[2] , y[1]),[1 2]) 
    end

    Q = (
        [w for w in W] # 入力データ
        .|> x-> Qlamda((V,x))  # Vに含まれるxのインデックスを探す 1が r=3,i=1にあるなら[3,1]が返る
        |> x -> tuple.(eachcol(x)...) # 上で帰ってくるのはMatrixなので行ごとにタプル化
        |> x -> filter(x->(N[x[1]] >= x[2]),x) # 余分な末尾の切り落とし
        |>Set  ) # 最後は集合にしておわり
    for (r,i) ∈  Q[1]
        println("r=",r,", i=",i)
    end

    
    M = 5000
    m = Model(HiGHS.Optimizer)
    
    set_attribute(m, "presolve", "on")
    set_attribute(m, "time_limit", 600.0)

    @variable(m, T, Int)
    @variable(m, s[r = 1:R, i = 1:N[r]] .>= 0, Int)
    @variable(m, z[r = 1:R, i = 1:N[r], k = 1:Capacity[V[r,i]]], Bin)
    @variable(m, y[r = 1:R, i = 1:N[r], u = 1:R, j = 1:N[u]], Bin)
    
    @variable(m, x[r = 1:R, i = 1:N[r], k = num_of_hoist], Bin)
    @constraint(m, 
                [r = 1:R, i = 1:N[r]],
                sum(x[r,i,k] for k in num_of_hoist) == 1)

    # #(5.3)
    # @constraint(m,
    #             [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u]; (r != u || i != j)], 
    #             s[u,j] - s[r,i] >= D[r,i] + E[V[r,i+1],V[u,j]] - M*(1-y[r,i,u,j]))
    
    # 途中経過
    # @constraint(m,
    #             [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u]; (r != u || i != j)], 
    #             s[u,j] - s[r,i] >= D[r,i] + E[V[r,i+1],V[u,j]] - M*(2-y[r,i,u,j] - x[r,i,1]))
    
    @constraint(m,
                [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u], k = num_of_hoist; (r != u || i != j)], 
                s[u,j] - s[r,i] >= D[r,i] + E[V[r,i+1],V[u,j]] - M*( 3 - y[r,i,u,j] - x[r,i,k] - x[u,j,k]))
    
    δ = 5 # r,iとu,jの搬送開始をどれだけずらすかのパラメータ
    @constraint(m,
                [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u], k = num_of_hoist; (r != u || i != j)], 
                s[u,j] - s[r,i] >= δ - M*( 3 - y[r,i,u,j] - x[r,i,k] - sum(x[u,j,l] for l in num_of_hoist if l != k)))
    
    #(5.4)
    @constraint(m, 
                [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u]; (r != u || i != j)],
                y[r,i,u,j] + y[u,j,r,i] == 1)
    
    #(5.5)
    @constraint(m, s[1,1] == 0)
    
    #(5.7)
    @constraint(m,
                [r =1:R, i = 1:N[r]],
                T >= s[r,i] + D[r,i] + E[V[r,i+1],1])

    #(5.8)
    @constraint(m,
                [r = 1:R, i = 2:N[r]],
                sum(z[r,i,k] for k = 1:Capacity[V[r,i]]) == 1)
    
    # (5.15)
    # L[r,i] <= s[r,i] - (s[r,i-1]+D[r,i-1]) + { ( ∑ k*z[r,i,k] ) - y[r,i-1,r,i] + 1 } *T <= U[r,i]
    # 線形化を考える
    # Φ = s[r,i] - (s[r,i-1]+D[r,i-1]) + { ( ∑ k*z[r,i,k] ) - y[r,i-1,r,i] + 1 } *T とすると
    # L[r,i] <= Φ <= U[r,i] なのでLとUで2つの式に分ける.
    # Φ の中で { ( ∑ k*z[r,i,k] ) - y[r,i-1,r,i] + 1 } *T が線形ではない.
    # 決定変数 T が決定変数 z と y に掛け算されてるので2次であり線形ではないため.
    # これの線形化のテクニックは別途調べる.

    # Σ の挙動を考えると Σⁿₖ₌₁  k=nの時も実行されるので、 for k = 1:0 にしてはいけない。 k = 1:0 だと実行されないため。
    # @constraint(m,
    # [r = 1:R, i = 2:N[r]],
    # sum([(k-1)* z[r,i,k] for k in 1:Capacity[V[r,i]]] ) - y[r,i,r,i-1]  >= 0)

    #(5.16) 
    @constraint(m,
                [r = 1:R, i = 2:N[r], k = 1:Capacity[V[r,i]]],
                s[r,i] - s[r,i-1] - D[r,i-1] + (k-1)*T >= L[r,i] - M*(2 - y[r,i-1,r,i] - z[r,i,k]))

    #(5.17)
    @constraint(m,
                [r = 1:R, i = 2:N[r], k = 1:Capacity[V[r,i]]],
                s[r,i] - s[r,i-1] - D[r,i-1] + (k-1)*T <= U[r,i] + M*(2 - y[r,i-1,r,i] - z[r,i,k]))
    #(5.18)
    @constraint(m,
                [r = 1:R, i = 2:N[r], k = 1:Capacity[V[r,i]]],
                s[r,i] - s[r,i-1] - D[r,i-1] + k*T >= L[r,i] - M*(1 + y[r,i-1,r,i] - z[r,i,k]))
    #(5.19)
    @constraint(m,
                [r = 1:R, i = 2:N[r], k = 1:Capacity[V[r,i]]],
                s[r,i] - s[r,i-1] - D[r,i-1] + k*T <= U[r,i] + M*(1 + y[r,i-1,r,i] - z[r,i,k]))

    #(5.22)
    @constraint(m,
                [w in W],
                sum(sum((k-1)*z[r,i,k] for k = 1:Capacity[w]) - y[r,i-1,r,i] + 1 for (r,i) in Q[w] if i > 1) <= Capacity[w])
    #(5.26)
    @constraint(m,
                [r = 1:R, i = 2:N[r]],
                sum(sum((k-1)*z[u,j,k] for k = 1:Capacity[V[r,i]]) + y[u,j,u,j-1] - y[u,j,r,i-1] for (u,j) in Q[V[r,i]] if j > 1) 
                + sum(y[u,j-1,r,i-1] for (u,j) in setdiff(Q[V[r,i]], Set(((r,i),))) if j > 1 ) 
                <= Capacity[V[r,i]]-1)
    #(5.27)
    @constraint(m,
                [r = 2:R, u = 2:R, i = 2:N[r], j = 2:N[u]; ((r != u || i != j) && V[r,i] == V[u,j])],
                y[r,i-1,r,i] + y[u,j-1,u,j] + y[r,i,u,j-1] + y[u,j,r,i-1] == 3)
    #(5.28)
    @constraint(m,
                [r = 1:R, i = 1:N[r]],
                s[r,i] >= 0)

    #(5.31)
    @objective(m, Min, T)
    JuMP.write_to_file(m, "model" * ".mps")
    optimize!(m)

    l(r,i) = sum([(k-1)*value(z[r,i,k]) for k in 1:Capacity[V[r,i]]] ) + value(y[r,i,r,i-1])
    t(r,i) = value(s[r,i]) - value(s[r,i-1]) - D[r,i-1] + l(r,i) * value(T)

    for ((r,i),) in s.data
        if i != 1
        println("l($(r),$(i)) = ", l(r,i))
        end
    end

    start_time = Dict([((r,i) => round(Int, value(s[r,i]))) for ((r,i),) in s.data])
    routes     = Dict([((r,i) => V[r,i]) for ((r,i),) in s.data])
    actual_time =  Dict([((r,i) =>  (i != 1) ? round(Int, t(r,i)) : 0) for ((r,i),) in s.data])

    y_values = Dict([((r,i,u,j) => value(y[r,i,u,j])) for ((r,i,u,j),) in y.data])

    for ((r,i),) in s.data
        if round(Int,value(x[r,i,1])) == 1
            println("transfer_resource[$(r),$(i),1] = ", value(transfer_resource[r,i,1]))
        end
    end
    transfer_no = Dict([((r,i,k) => round(Int,value(x[r,i,k]))) for ((r,i,k),) in x.data])
    # println(value.(x))

    # sortを辞書型に変換して返してはいけない。順序が変わってしまうため。
    return sort(start_time; byvalue = true), routes, actual_time, y_values, value(T), transfer_no
end

struct cyclic_hoist_schedule
    hoist_move_starttime::Dict{Tuple{Int,Int},Int}
    to::Dict{Tuple{Int,Int},Int}
    from::Dict{Tuple{Int,Int},Int}
    has_stuff::Dict{Tuple{Int,Int},Bool}
end

function main(input::cyclic_hoist_schedule) 
    main(W,R,N,Capacity,V,L,U,E,D)
end

"""
    load_input_data(filepath::String) → params_dict

JSONファイルから入力データを読み込み、辞書にまとめたものを返す。
"""
function load_input_data(filepath::String)
    raw = JSON.parsefile(filepath)

    W = Int.(raw["W"])
    R = raw["R"]
    N = Int.(raw["N"])

    Capacity = ones(Int, maximum(W))
    for (k, v) in raw["Capacity"]
        Capacity[parse(Int, k)] = v isa Integer ? Int(v) : parse(Int, string(v))
    end

    V = Int.(hcat(raw["V"]...))
    L = Int.(hcat(raw["L"]...))
    U = Int.(hcat(raw["U"]...))

    max_w = maximum(W)
    E = fill(2, (max_w + 1, max_w + 1))
    for i in W, j in i:max_w+1
        E[j, j] = 0
        if j > i
            E[i, j] = sum(E[k, k+1] for k in i:j-1)
            E[j, i] = E[i, j]
        end
    end

    D = zeros(Int, (R, max_w + 1))
    for r in 1:R, i in 1:N[r]
        D[r, i] = E[V[r, i], V[r, i+1]] + 20
        if ( V[r, i] == V[r, i+1])
            D[r, i] = 0
        end
    end

    return W, R, N, Capacity, V, L, U, E, D
end

"""
    read_json(filepath::String)

JSONファイルからデータを読み込んでmain関数を呼び出す。
"""
function read_json(filepath::String)
    W, R, N, Capacity, V, L, U, E, D = load_input_data(filepath)
    return W, R, N, Capacity, V, L, U, E, D
end

#############################################################
# デフォルトのデータ定義
#############################################################
# コマンドライン引数で入力ファイルを指定する
if length(ARGS) < 1
    println("使用方法: julia hoist_cyclic.jl <input_file.json>")
    input_file = "input_data.json"  # デフォルトの入力ファイル名
    # exit(1)
else
    input_file = ARGS[1]
end
println("入力ファイル: ", input_file)

W, R, N, Capacity, V, L, U, E, D = read_json(input_file)
st, rt, at, y_val, cycletime, transfer_resource = main(W,R,N,Capacity,V,L,U,E,D)

#############################################################

"""
    plot_schedule(st, rt, at, V, L, U, E, D, N, cycletime; colors = [:red, :blue, :green])

スケジューリング結果をプロットする関数。

# Argumentst
- `st`: 開始時間のDictまたはOrderedDict
- `rt`: ルートのDict
- `at`: 実時間Dict
- `V`: 作業スケジュール行列
- `L`: 下限時間制約
- `U`: 上限時間制約
- `E`: 移動時間行列
- `D`: 作業時間行列
- `N`: 各ロボットの仕事数
- `cycletime`: サイクルタイム
- `colors`: ロボットごとの色（オプション）

# Returns
- `fig`: MakieのFigureオブジェクト
"""
function plot_schedule(
        st,
        rt::Dict{Tuple{Int,Int},Int},
        at,
        V::Matrix{Int},
        L::Matrix{Int},
        U::Matrix{Int},
        E::Matrix{Int},
        D::Matrix{Int},
        N::Vector{Int},
        cycletime,
        transfer_resource;
        colors = [:red, :blue, :green]
    )

    fig = Figure()

    ax = Axis(fig[1,1])
    ax2 = Axis(fig[2,1])
    ax3 = Axis(fig[3,1])
    axbottom = Axis(fig[4, 1], yticks = ([1,2], ["",""]), ylabel = "Handling")
    linkxaxes!(ax, axbottom)
    linkxaxes!(ax, ax2)
    linkxaxes!(ax, ax3)

    rz1 = 1
    iz1 = 1

    for (count, (r, i)) in enumerate(keys(st))
        if transfer_resource[r,i,1] == 1
            if (i <= N[r])
                lines!(ax, [st[r,i], st[r,i] + D[r,i]], [V[r,i], V[r,i+1]]; color = :black)
                text!(ax, st[r,i] + D[r,i], V[r,i+1], text="$(st[r,i] + D[r,i])", align = (:left, :bottom))
                
                barplot!(axbottom, 1, st[r,i] + D[r,i], fillto = st[r,i], direction = :x, color = colors[r])
            end

            if ((r == rz1) && abs(i - iz1) > 1 && iz1 <= N[rz1]) || (r != rz1) && iz1 <= N[rz1] && count != 1 && count != length(st)
                lines!(ax, [st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], st[r,i]], [V[rz1,iz1+1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 1, st[r,i], fillto = st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], direction = :x, color = colors[rz1])
            elseif iz1 == N[rz1]
                lines!(ax, [st[rz1,iz1], st[r,i]], [V[rz1,iz1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 1, st[r,i], fillto = st[rz1,iz1], direction = :x, color = colors[rz1])
            end

            if (i <= N[r])
                println("actual_time[$(r),$(i)] = ", at[r,i])
                if cycletime > (st[r,i] - at[r,i]) && (st[r,i] - at[r,i]) > 0
                    lines!(ax, [st[r,i] - at[r,i], st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                else
                    lines!(ax, [0, st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                    lines!(ax, [cycletime + (st[r,i] - at[r,i]), cycletime], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                end
            end

            text!(ax, st[r,i], V[r,i], text="r$(r),i$(i) $(st[r,i])", align = (:left, :top))

            rz1 = r
            iz1 = i
        end
    end

    rz1 = 1
    iz1 = 1

if (length(keys(transfer_resource)) != length(keys(st))) 
    for (count, (r, i)) in enumerate(keys(st))
        if transfer_resource[r,i,2] == 1
            if (i <= N[r])
                lines!(ax2, [st[r,i], st[r,i] + D[r,i]], [V[r,i], V[r,i+1]]; color = :black)
                text!(ax2, st[r,i] + D[r,i], V[r,i+1], text="$(st[r,i] + D[r,i])", align = (:left, :bottom))
                
                barplot!(axbottom, 2, st[r,i] + D[r,i], fillto = st[r,i], direction = :x, color = colors[r])
            end

            if ((r == rz1) && abs(i - iz1) > 1 && iz1 <= N[rz1]) || (r != rz1) && iz1 <= N[rz1] && count != 1 && count != length(st)
                lines!(ax2, [st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], st[r,i]], [V[rz1,iz1+1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 2, st[r,i], fillto = st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], direction = :x, color = colors[rz1])
            elseif iz1 == N[rz1]
                lines!(ax2, [st[rz1,iz1], st[r,i]], [V[rz1,iz1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 2, st[r,i], fillto = st[rz1,iz1], direction = :x, color = colors[rz1])
            end

            if (i <= N[r])
                println("actual_time[$(r),$(i)] = ", at[r,i])
                if cycletime > (st[r,i] - at[r,i]) && (st[r,i] - at[r,i]) > 0
                    lines!(ax2, [st[r,i] - at[r,i], st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                else
                    lines!(ax2, [0, st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                    lines!(ax2, [cycletime + (st[r,i] - at[r,i]), cycletime], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                end
            end

            text!(ax2, st[r,i], V[r,i], text="r$(r),i$(i) $(st[r,i])", align = (:left, :top))

            rz1 = r
            iz1 = i
        end
    end
end
#############################################################################################################
    rz1 = 1
    iz1 = 1
    for (count, (r, i, k)) in enumerate(keys(sort(transfer_resource; byvalue = true)))
        if k == 1
            if (i <= N[r])
                lines!(ax3, [st[r,i], st[r,i] + D[r,i]], [V[r,i], V[r,i+1]]; color = :black)
                text!(ax3, st[r,i] + D[r,i], V[r,i+1], text="$(st[r,i] + D[r,i])", align = (:left, :bottom))
                
                barplot!(axbottom, 1, st[r,i] + D[r,i], fillto = st[r,i], direction = :x, color = colors[r])
            end

            if ((r == rz1) && abs(i - iz1) > 1 && iz1 <= N[rz1]) || (r != rz1) && iz1 <= N[rz1] && count != 1 && count != length(st)
                lines!(ax3, [st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], st[r,i]], [V[rz1,iz1+1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 1, st[r,i], fillto = st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], direction = :x, color = colors[rz1])
            elseif iz1 == N[rz1]
                lines!(ax3, [st[rz1,iz1], st[r,i]], [V[rz1,iz1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 1, st[r,i], fillto = st[rz1,iz1], direction = :x, color = colors[rz1])
            end

            if (i <= N[r])
                println("actual_time[$(r),$(i)] = ", at[r,i])
                if cycletime > (st[r,i] - at[r,i]) && (st[r,i] - at[r,i]) > 0
                    lines!(ax3, [st[r,i] - at[r,i], st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                else
                    lines!(ax3, [0, st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                    lines!(ax3, [cycletime + (st[r,i] - at[r,i]), cycletime], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                end
            end

            text!(ax3, st[r,i], V[r,i], text="r$(r),i$(i) $(st[r,i])", align = (:left, :top))

            rz1 = r
            iz1 = i
        end
    end

    rz1 = 1
    iz1 = 1

    for (count, (r, i, k)) in enumerate(keys(sort(transfer_resource; byvalue = true)))
        if k == 2
            if (i <= N[r])
                lines!(ax3, [st[r,i], st[r,i] + D[r,i]], [V[r,i], V[r,i+1]]; color = :black)
                text!(ax3, st[r,i] + D[r,i], V[r,i+1], text="$(st[r,i] + D[r,i])", align = (:left, :bottom))
                
                barplot!(axbottom, 2, st[r,i] + D[r,i], fillto = st[r,i], direction = :x, color = colors[r])
            end

            if ((r == rz1) && abs(i - iz1) > 1 && iz1 <= N[rz1]) || (r != rz1) && iz1 <= N[rz1] && count != 1 && count != length(st)
                lines!(ax3, [st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], st[r,i]], [V[rz1,iz1+1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 2, st[r,i], fillto = st[rz1,iz1]+E[V[rz1,iz1+1],V[rz1,iz1+1]]+D[rz1,iz1], direction = :x, color = colors[rz1])
            elseif iz1 == N[rz1]
                lines!(ax3, [st[rz1,iz1], st[r,i]], [V[rz1,iz1], V[r,i]], linestyle = :dash; color = :black)
                barplot!(axbottom, 2, st[r,i], fillto = st[rz1,iz1], direction = :x, color = colors[rz1])
            end

            if (i <= N[r])
                println("actual_time[$(r),$(i)] = ", at[r,i])
                if cycletime > (st[r,i] - at[r,i]) && (st[r,i] - at[r,i]) > 0
                    lines!(ax3, [st[r,i] - at[r,i], st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                else
                    lines!(ax3, [0, st[r,i] ], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                    lines!(ax3, [cycletime + (st[r,i] - at[r,i]), cycletime], [V[r,i] + (r-1)*0.025, V[r,i] + (r-1)*0.025]; color = colors[r]) 
                end
            end

            text!(ax3, st[r,i], V[r,i], text="r$(r),i$(i) $(st[r,i])", align = (:left, :top))

            rz1 = r
            iz1 = i
        end
    end
##################################################################################################

    rowsize!(fig.layout, 4, Auto(0.1))
    
    return fig
end

# 描画処理を関数として呼び出す
fig = plot_schedule(st, rt, at, V, L, U, E, D, N, cycletime, transfer_resource)
