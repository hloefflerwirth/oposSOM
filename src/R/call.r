pipeline.call <- function(fn, env)
{
  environment(fn) <- env
  return(fn())
}
