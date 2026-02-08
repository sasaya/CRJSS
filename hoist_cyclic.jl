using JuMP, HiGHS, GLMakie #CairoMakie


struct resultHoist
    hoist_move_starttime
    to
    from
    has_stuff
end

function main(W,R,N,Capacity,V,L,U,E,D)
    # W = 12
    # R = 3
    # N = [12, 8, 8]

    # Capacity = ones(Int,W)
    # Capacity[2] = 2
    # Capacity[7] = 3
    # Capacity[8] = 2

    # V = [ 1   2   3   4   5   6   7    8    9  10  11  12;
    #       1   2   4   5   7   8   11   12   0  0   0   0 ;
    #       1   3   5   7   8   9   10   12   0  0   0   0 ]

    # L = [ 30 150  60  60  30  40 350  90  70  45  60  30;
    #       30 130  30  40 420  50  90  70   0   0   0   0;
    #       30  80  40 300  50  90  90  70   0   0   0   0]

    # U = [ 60 350  90 120  75 120 800 160 200  90 120  80;
    #       50 280  70  90 850 120 160 200   0   0   0   0;
    #       50 220  90 550 120 150 160 200   0   0   0   0]

    function Qlamda(y) 
        return getindex.(findall(x -> x == y[2] , y[1]),[1 2]) 
    end

    Q = [Qlamda((V,x)) |> x -> tuple.(eachcol(x)...) |> Set  for x in W]  
    for (r,i) ∈  Q[1]
        println("r=",r,", i=",i)
    end
     
    # E = maximum(W)+1 |> x -> zeros(Int,(x,x))
    # for i in W, j in i:maximum(W)+1
    #     E[i,i+1] = 2
    #     if (j > i)
    #         E[i,j] = 2 * ((j-1) - i)
    #         E[j,i] = E[i,j]
    #     end
    # end

    # D = zeros(Int,(R,maximum(W)+1))
    # for r in 1:R, i in 1:N[r]-1
    #     D[r,i] = E[V[r,i], V[r,i+1]] + 20
    # end

    M = 5000
    m = Model(HiGHS.Optimizer)
    
    set_attribute(m, "presolve", "on")
    set_attribute(m, "time_limit", 36000.0)

    @variable(m, T, Int)
    @variable(m, s[r = 1:R, i = 1:N[r]] .>= 0, Int)
    # @variable(m, t[r = 1:R, i = 1:N[r]], Int)
    @variable(m, z[r = 1:R, i = 1:N[r], k = 1:Capacity[V[r,i]]], Bin)
    @variable(m, y[r = 1:R, i = 1:N[r], u = 1:R, j = 1:N[u]], Bin)
    
    
    #(5.3)
    @constraint(m,
                [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u]; (r != u || i != j)], 
                s[u,j] - s[r,i] >= D[r,i] + E[V[r,i+1],V[u,j]] - M*(1-y[r,i,u,j]))
    #(5.4)
    @constraint(m, 
                [r = 1:R, u = 1:R, i = 1:N[r], j = 1:N[u]; (r != u || i != j)],
                y[r,i,u,j] + y[u,j,r,i] == 1)
    
    #(5.5)
    @constraint(m, s[1,1] == 0)
    
    #(5.6)
    # @constraint(m, s[r,i] >= D[1,1])

    #(5.7)
    @constraint(m,
                [r =1:R, i = 1:N[r]],
                T >= s[r,i] + D[r,i] + E[V[r,i+1],1])

    #(5.8) zのindexとしてのkと, サイクルTの係数としてのkの値は1つずれる(Juliaのindexが1始まりだから.)
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
    #(5.29)

    #(5.30)

    #(5.31)
    @objective(m, Min, T)
    
    optimize!(m)

    # start_time = Array{Tuple{Int64, Int64}, Int64}(R,W)
                     # Dcit(zip(a,b))でもいい

    l(r,i) = sum([(k-1)*value(z[r,i,k]) for k in 1:Capacity[V[r,i]]] ) + value(y[r,i,r,i-1])
    t(r,i) = value(s[r,i]) - value(s[r,i-1]) - D[r,i-1] + l(r,i) * value(T)

    for ((r,i),) in s.data
        if i != 1
        println("l($(r),$(i)) = ", l(r,i))
        end
    end

    start_time = Dict([((r,i) => round(value(s[r,i]))) for ((r,i),) in s.data])
    routes     = Dict([((r,i) => V[r,i]) for ((r,i),) in s.data])
    actual_time =  Dict([((r,i) =>  (i != 1) ? round(t(r,i)) : 0) for ((r,i),) in s.data])

    y_values = Dict([((r,i,u,j) => value(y[r,i,u,j])) for ((r,i,u,j),) in y.data])
    
    # sort(start_time; byvalue = true)
    return sort(start_time; byvalue = true), routes, actual_time, y_values, value(T)
end


#############################################################
W = 1:12
R = 3
N = [12, 8, 8]

Capacity = ones(Int,W)
Capacity[2] = 2
Capacity[7] = 3
Capacity[8] = 2

# V = [ 1   2   3   4   5   6   7    8    9  10  11  12;
#       1   2   4   5   7   8   11   12   0  0   0   0 ;
#       1   3   5   7   8   9   10   12   0  0   0   0 ]

V = [1   2   3   4   5   6   7    8    9  10  11  12 13;
     1   2   4   5   7   8   11   12   13  0   0   0  0;
     1   3   5   7   8   9   10   12   13  0   0   0  0]

L = [0 30 150  60  60  30  40 350  90  70  45  60  30  0;
     0 30 130  30  40 420  50  90  70   0   0   0   0  0;
     0 30  80  40 300  50  90  90  70   0   0   0   0  0]

U = [0 60 350  90 120  75 120 800 160 200  90 120  80  0;
     0 50 280  70  90 850 120 160 200   0   0   0   0  0;
     0 50 220  90 550 120 150 160 200   0   0   0   0  0]

st, rt, at, y_val, cycletime = main(W,R,N,Capacity,V,L,U,E,D)

E = maximum(W)+1 |> x -> zeros(Int,(x,x))
for i in W, j in i:maximum(W)+1
    E[i,i+1] = 2
    if (j > i)
        E[i,j] = 2 * ((j-1) - i)
        E[j,i] = E[i,j]
    end
end

D = zeros(Int,(R,maximum(W)+1))
for r in 1:R, i in 1:N[r]-1
    D[r,i] = E[V[r,i], V[r,i+1]] + 20
end
#############################################################

fig = Figure()

ax = Axis(fig[1,1])

# lines!(ax,
#         get.(Ref(st), keys(st), missing),
#         get.(Ref(rt), keys(st), missing)
#     )
rz1 = 1
iz1 = 1

c = [:red, :blue, :green]

for (count,(r,i)) in enumerate(keys(st))
    if (i != N[r])
    # 同じパーツが次のところに移動するときは実線
        lines!(ax, [st[r,i], st[r,i] + D[r,i]], [rt[r,i], rt[r,i+1]]; color = :black)
        text!(ax,st[r,i] + D[r,i], rt[r,i+1], text="$(st[r,i] + D[r,i])", align = (:left, :bottom))
    end

    # 空荷移動は破線, これ以降は時間でソートされてる前提. unloadした後に別のものをloadするときを想定しているため.
    if  ((r == rz1) && abs(i - iz1) > 1 && iz1 != N[rz1]) || (r != rz1) && iz1 != N[rz1] && count != 1 && count != length(st)
        lines!(ax, [st[rz1,iz1]+E[rt[rz1,iz1+1],rt[rz1,iz1+1]]+D[rz1,iz1], st[r,i]], [rt[rz1,iz1+1], rt[r,i]], linestyle = :dash; color = :black)
    elseif i == N[r]
        lines!(ax, [st[r,i], st[r,i]+E[V[r,i],V[r,1]]], [rt[r,i], rt[r,1]], linestyle = :dash; color = :black)
    end

    # 横棒は色付き
    if (i != N[r])
        println("actual_time[$(r),$(i)] = ", at[r,i])
    # 同じパーツが次のところに移動するときは実線
    if cycletime > (st[r,i] - at[r,i]) && (st[r,i] - at[r,i]) > 0
        lines!(ax, [st[r,i] - at[r,i], st[r,i] ], [rt[r,i] + (r-1)*0.025, rt[r,i] + (r-1)*0.025]; color = c[r]) 
    else
        lines!(ax, [0, st[r,i] ], [rt[r,i] + (r-1)*0.025, rt[r,i] + (r-1)*0.025]; color = c[r]) 
        lines!(ax, [cycletime + (st[r,i] - at[r,i]), cycletime], [rt[r,i] + (r-1)*0.025, rt[r,i] + (r-1)*0.025]; color = c[r]) 
    end
    end

    text!(ax,st[r,i], rt[r,i], text="r$(r),i$(i) $(st[r,i])", align = (:left, :bottom))
    global rz1 = r
    global iz1 = i

end
fig
