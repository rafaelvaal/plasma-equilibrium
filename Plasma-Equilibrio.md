```julia
#Valores en la frontera:
α = 1.0   #valor en r = 0.01
β = 0.0 #valor de y en la superficie del cilíndro
a = 0.01 #frontera izquierda, centro del cilíndro. No se toma exactamente en el cero para evitar a la singularidad
b = 1.0 #frontera exterior del cilindro
M = 100 #número máximo de iteraciones
tolerancia = 0.01
N = 100 #número de subintervalos
h = (b-a)/N
##################################
k = 1
tk = (β-α)/(b-a) #aproximación inicial de la pendiente a variar

#Definición de Parámetros:
A4 = -0.18
A3 = -0.0
A2 = 0.0
A1 = -0.0
A0 = 7
a0 =0/2
a1 =A1/2
a2 =A2/2
a3 =A3/2
a4 =A4/2
b0 = A0
b1 = A1/2
b2 = A2/2
b3 = A3/2
b4 = A4/2

sol = [(a,α)]                  ##Éstas listas guardarán las coordenadas
pres = [(a,-(a4*α^4+a3*α^3+a2*α^2+a1*α+a0))]
paramb = [(a,-2*sqrt((b4*α^4+b3*α^3+b2*α^2+b1*α+b0)))]

function f(x,y,y2)
    global A1, A2, A3, A4, A0
    f = -(1/x)*y2 + A4*y^4+A3*y^3+A2*y^2+A1*y+A0
    return f
end
function fy(x,y,y2)
    global A1, A2, A3, A4, A0
    f = 4*A4*y^3 + 3*A3*y^2 + 2*A2*y + A1
    return f
end
function fy2(x,y,y2)########derivada de f respecto a y'
    f = -(1/x)
    return f
end

while k<=M
global α,tk,u1,u2
w1 = [α]
w2 = [tk]
u1 = 0
u2 = 1
    for i in 2:N
        global u1,u2
        x = a + (i-2)*h
##### Definición de términos de Runge-Kutta
        k11 = h*w2[i-1]
        k12 = h*f(x,w1[i-1],w2[i-1])
        k21 = h*(w2[i-1]+(0.5)*k12)
        k22 = h*f(x+(0.5)*h, w1[i-1]+(0.5)*k11, w2[i-1]+(0.5)*k12)
        k31 = h*(w2[i-1]+(0.5)*k22)
        k32 = h*f(x+(0.5)*h, w1[i-1]+(0.5)*k21, w2[i-1]+(0.5)*k22)
        k41 = h*(w2[i-1]+(0.5)*k32)
        k42 = h*f(x+h, w1[i-1]+(0.5)*k31, w2[i-1]+(0.5)*k32)

        w1y = w1[i-1] + (k11 + 2*k21 + 2*k31 + k41)*(1/6)
        w2y = w2[i-1] + (k12 + 2*k22 + 2*k32 + k42)*(1/6)
        push!(w1, w1y)
        push!(w2, w2y)
        c11 = h*u2
        c12 = h*(fy(x,w1[i-1],w2[i-1])*u1 + fy2(x,w1[i-1],w2[i-1])*u2)
        c21 = h*(u2+(0.5)*c12)
        c22 = h*(fy(x+(0.5)*h, w1[i-1], w2[i-1])*(u1+(0.5)*c11) + fy2(x+(0.5)*h, w1[i-1], w2[i-1])*(u2+(0.5)*c21))
        c31 = h*(u2+(0.5)*c22)
        c32 = h*(fy(x+(0.5)*h, w1[i-1], w2[i-1])*(u1+(0.5)c21) + fy2(x+(0.5)*h, w1[i-1], w2[i-1])*(u2+(0.5)*c22) )
        c41 = h*(u2+c32)
        c42 = h*(fy(x+(0.5)*h, w1[i-1], w2[i-1])*(u1+c31) + fy2(x+(0.5)*h, w1[i-1], w2[i-1])*(u2+c32))
        u1 = u1 + (1/6)*(c11+2*c21+2*c31+c41)
        u2 = u2 + (1/6)*(c12+2*c22+2*c32+c42)
    end
        if abs(w1[N]-β)<=tolerancia
            for i in 1:N
                #global sol
                x = a + i*h
                push!(sol,(x,w1[i]))
                P = -(a4*w1[i]^4+a3*w1[i]^3+a2*w1[i]^2+a1*w1[i]+a0)
                push!(pres,(x,P))
                b = -2*sqrt((b4*w1[i]^4+b3*w1[i]^3+b2*w1[i]^2+b1*w1[i]+b0))
                push!(paramb,(x,b))
            end
            break
        else
            global u1,k,β
            tk = tk - (w1[N]-β)/u1
            k = k + 1
        end
end
#popfirst!(sol)
#popfirst!(pres)
#popfirst!(paramb)
#using Plots
#pyplot()
#plot(sol, lw=2, title="Flujo Magnético", xaxis=("Radio r"), yaxis=("Flujo Magnético ψ"), label = "Flujo ψ")
#plot!(pres,  lw=2,  xaxis=("Radio r"), yaxis=("μP"),title="Perfil de Presión", label = "Presión P")
#plot(paramb, lw=2, xaxis=("Radio r"), yaxis=("-1/2b^2"),title="Perfiles de Flujo, Parámetro b y Campo Magnético", label = "-1/2b^2")
#plot!(paramb, lw=2, xaxis=("Radio r"), yaxis=("b"),title="Soluciones al Parámetro b", legend = false)
#savefig("paramb")

```


```julia

```
