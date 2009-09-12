`error` <- function(forecast, forecastbench, true, method = c("me", "mae", "mse", "sse", "rmse", "mdae", "mdse", "mape", "mdape", 
            "smape", "smdape", "rmspe", "rmdspe", "mrae", "mdrae", "gmrae", "relmae", "relmse", "mase", "mdase", "rmsse"))
{
  if (method == "me"){
      val = me(forecast,true)
  }
  if (method == "mae"){
      val = mae(forecast,true)
  }
  if (method == "mse"){
      val = mse(forecast,true)
  }
  if (method == "sse"){
      val = sse(forecast,true)
  }
  if (method == "rmse"){
      val = rmse(forecast,true)
  }
  if (method == "mdae"){
      val = mdae(forecast,true)
  }
  if (method == "mdse"){
      val = mdse(forecast,true)
  }
  if (method == "mape"){
      val = mape(forecast,true)
  }
  if (method == "mdape"){
      val = mdape(forecast,true)
  }
  if (method == "smape"){
      val = smape(forecast,true)
  }
  if (method == "smdape"){
      val = smdape(forecast,true)
  }
  if (method == "rmspe"){
      val = rmspe(forecast,true)
  }
  if (method == "rmdspe"){
      val = rmdspe(forecast,true)
  }
  if (method == "mrae"){
      val = mrae(forecast,forecastbench,true)
  }
  if (method == "mdrae"){
      val = mdrae(forecast,forecastbench,true)
  }
  if (method == "gmrae"){
      val = gmrae(forecast,forecastbench,true)
  }
  if (method == "relmae"){
      val = relmae(forecast,forecastbench,true)
  }
  if (method == "relmse"){
      val = relmse(forecast,forecastbench,true)
  }
  if (method == "mase"){
      val = mase(forecast,true)
  }
  if (method == "mdase"){
      val = mdase(forecast,true)
  }
  if (method == "rmsse"){
      val = rmsse(forecast,true)
  }
  return(val)
}


