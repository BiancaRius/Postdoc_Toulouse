library(rcontroll)

#load data
data("TROLLv3_species")
data("TROLLv3_climatedaytime12")
data("TROLLv3_daytimevar")
data("TROLLv3_input")

#run the simulation
troll(
  name = "test_outputs",
  path = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/rconTROLL/",
  global = generate_parameters(
    cols = 100, rows = 100,
    iterperyear = 12, nbiter = 1200
  ),
  species = TROLLv3_species,
  climate = TROLLv3_climatedaytime12,
  daily = TROLLv3_input
)


#load the outputs and create a trollsim() object
t1 = load_output("test_outputs","/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/rconTROLL/test_outputs/")

#The plot is supper weird
autoplot(t1)
