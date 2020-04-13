options(digits = 4)
solve_x <- function(x0, x1, y0, y1) {
  return (x1 - (y1 * (x1 - x0)) / (y1 - y0))
}
solve_ea <- function(prev_x, x) {
  return (abs((x - prev_x) / x) * 100)
}

SecantMethod <- function(f,
                         given_x0,
                         given_x1,
                         macheps,
                         max,
                         verbose) {
  if (missing(verbose)) {
    verbose = TRUE
  }
  
  x0 = given_x0
  x1 = given_x1
  y0 = NA
  y1 = NA
  x = NA
  y = NA
  prev_x = NA
  rownames = c(1:max)
  colnames = c("x0", "x1", "y0", "y1", "x", "y", "ea")
  m = matrix(
    dimnames = list(rownames, colnames),
    nrow = max,
    ncol = 7
  )
  
  num_of_iter = NA
  for (i in 1:max) {
    num_of_iter = i
    y0 = f(x0)
    y1 = f(x1)
    x = solve_x(x0, x1, y0, y1)
    y = f(x)
    ea = solve_ea(prev_x, x)
    m[i, ] = c(x0, x1, y0, y1, x, y, ea)
    x0 = x1
    x1 = x
    prev_x = x
    if (y==0){
      break
    }
    if (!is.na(ea) && ea < macheps) {
      break
    }
  }
  if (verbose) {
    print(m[1:num_of_iter,])
  }
  return(list(x = m[num_of_iter, "x"], iterations = num_of_iter,ea = ea))
}

MullerMethod <- function(f,
                         given_x0,
                         given_x1,
                         given_x2,
                         macheps,
                         max,
                         verbose) {
  if (missing(verbose)) {
    verbose = TRUE
  }
  rownames = c(1:max)
  colnames = c("x0", "x1", "x2", "y0", "y1", "y2", "A", "B", "C", "x3", "y3", "ea")
  m = matrix(
    dimnames = list(rownames, colnames),
    nrow = max,
    ncol = 12
  )
  x0 = given_x0
  x1 = given_x1
  x2 = given_x2
  
  num_of_iter = NA
  
  for (i in 1:max) {
    num_of_iter = i
    y0 = f(x0)
    y1 = f(x1)
    y2 = f(x2)
    h0 = x1 - x0
    h1 = x2 - x1
    d0 = (y1 - y0) / h0
    d1 = (y2 - y1) / h1
    A = (d1 - d0) / (h1 + h0)
    B = A * h1 + d1
    C = y2
    
    positive_form =  B + sqrt(B ** 2 - 4 * A * C)
    negative_form = B - sqrt(B ** 2 - 4 * A * C)
    
    sign = if (abs(negative_form) < abs(positive_form))
      positive_form
    else
      negative_form
    x3 = x2 - 2 * C / sign

    y3 = f(x3)

    ea = solve_ea(x2, x3)

    m[i,] = c(x0, x1, x2, y0, y1, y2, A, B, C, x3, y3, ea)
    x0 = x1
    x1 = x2
    x2 = x3
    if (ea < macheps) {
      print(ea<macheps)
      break
    }
  }

  if (verbose) {
    print(m[1:num_of_iter,])
  }
  return(list(
      given_x0 = given_x0,
      given_x1 = given_x1,
      given_x2 = given_x2,
      x3 = x3,
      iterations = num_of_iter,
      ea = ea
    ))
}

f <- function(x) {
  return (sin(x) + cos(1 + x**2)-1)
}

s = SecantMethod(f,1,3,1*(10**-9),100000,FALSE)

ff <- function(x){
  return(x**3+3.5*x**2-40)
}
m = MullerMethod(ff, 1, 2.5, 4, 1 * (10 ** -9), 100000)