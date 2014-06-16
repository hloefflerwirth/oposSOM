
Get.Area = function(x, y)
{
  o = order(x)

  x = x[o]
  y = y[o]

  A = 0

  for (i in seq_len(length(x)-1))
  {
    dx = abs(x[i+1] - x[i])

    A = A + (dx * y[i])
  }

  A = A + (dx * y[length(x)])

  return(A)
}

