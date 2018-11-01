function euler(f, t₀, y₀, tₙ; n = 100)
  h = (tₙ - t₀) / n
  t = linspace(t₀, tₙ, n + 1)
  y = zeros(n + 1)
  y[1] = y₀
  for i = 1:n
    y[i+1] = y[i] + h * f(t[i], y[i])
  end
  return t, y
end

function euler_aperfeicoado(f,t_0,y_0,t_n;n=100)
    h = (t_n-t_0)/n
    t = linspace(t_0,t_n,n+1)
    y = zeros(n+1)
    y[1] = y_0
    alpha = 1/2
    beta = 1/2
    gamma = 1
    delta = 1
    for i=1:n
        k1 = f(t[i],y[i])
        k2 = f(t[i] + delta*h,y[i] + gamma*h*k1)
        y[i+1] = y[i]+h*(alpha*k1 + beta*k2)
    end  
    return t,y
end

function heun(f,t_0,y_0,t_n;n=100)
    h = (t_n-t_0)/n
    t = linspace(t_0,t_n,n+1)
    y = zeros(n+1)
    y[1] = y_0
    alpha = 1/4
    beta = 3/4
    gamma = 2/3
    delta = 2/3
    for i=1:n
        k1 = f(t[i],y[i])
        k2 = f(t[i] + delta*h,y[i] + gamma*h*k1)
        y[i+1] = y[i]+h*(alpha*k1 + beta*k2)
    end  
    return t,y
end

function midpoint(f,t_0,y_0,t_n;n=100)
    h = (t_n-t_0)/n
    t = linspace(t_0,t_n,n+1)
    y = zeros(n+1)
    y[1] = y_0
    alpha = 0
    beta = 1
    gamma = 1/2
    delta = 1/2
    for i=1:n
        k1 = f(t[i],y[i])
        k2 = f(t[i] + delta*h,y[i] + gamma*h*k1)
        y[i+1] = y[i]+h*(alpha*k1 + beta*k2)
    end  
    return t,y
end

function rungekutta4(f, t_0, y_0, t_n; n=100)
    h = (t_n-t_0)/n
    t = linspace(t_0,t_n,n+1)
    y = zeros(n+1)
    y[1] = y_0
    for i=1:n
        k1 = f(t[i], y[i])
        k2 = f(t[i]+h/2,y[i]+h*k1/2)
        k3 = f(t[i]+h/2,y[i]+h*k2/2)
        k4 = f(t[i]+h, y[i]+k3*h)
        y[i+1] = y[i]+(k1+2*k2+2*k3+k4)*h/6
    end
    return t, y    
end
