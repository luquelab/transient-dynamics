# 18-10-2023. This is a transcription from Toni's colab notebook: https://colab.research.google.com/github/luquelab/tmp/blob/main/dynamic_regimes.ipynb#scrollTo=rs6EhpOyXYKR

# Imports
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import json
#from google.colab import files


# Definitions

##############
## Class to generate objects containing the necessary information about a model
##############
class Model:
  def __init__(self, observer:dict={}, variables:dict={}, parameters:dict={}):
    self.observer = observer
    self.variables = variables
    self.parameters = parameters

  def __str__(self):
    return json.dumps(self.__dict__)

  def __repr__(self):
    return self.__str__()

  # Add variables
  def add_variables(self, variables:list):
    for variable in variables:
      self.variables[variable] = {}

  # Add symbols to variables
  def add_symbols_to_variables(self,variables:list, symbols:list):
    for index, variable in enumerate(variables):
      self.variables[variable]['symbol'] = symbols[index]

  # Add units to variables
  def add_units_to_variables(self, variables:list, units:list):
    for index, variable in enumerate(variables):
      self.variables[variable]['unit'] = units[index]

  # Add status to variables
  def add_status_to_variables(self,variables:list, independence:list):
    for index,variable in enumerate(variables):
      self.variables[variable]['independent'] = independence[index]

  # Add mechanisms contributing to the rate of change of variables
  def add_mechanisms(self, variables:list, mechanisms:list, expressions:list):
    for index,variable in enumerate(variables):
      self.variables[variable]['mechanisms'] = {}

      variable_mechanisms = mechanisms[index]
      variable_mechanisms_expressions = expressions[index]
      for index,mechanism in enumerate(variable_mechanisms):
         self.variables[variable]['mechanisms'][mechanism] = variable_mechanisms_expressions[index]

  # Add initial values to variables
  def add_initial_values_to_variables(self,variables:list,values:list):
    for index,variable in enumerate(variables):
      self.variables[variable]['initial_value'] = values[index]

  # Add parameters
  def add_parameters(self, parameters:list):
    for parameter in parameters:
      self.parameters[parameter] = {}

  # Add symbols to parameters
  def add_symbols_to_parameters(self,parameters:list, symbols:list):
    for index, parameter in enumerate(parameters):
      self.parameters[parameter]['symbol'] = symbols[index]

  # Add units to parameters
  def add_units_to_parameters(self, parameters:list, units:list):
    for index, parameter in enumerate(parameters):
      self.parameters[parameter]['unit'] = units[index]

  # Add initial values to parameters
  def add_initial_values_to_parameters(self,parameters:list,values:list):
    for index,parameter in enumerate(parameters):
      self.parameters[parameter]['value'] = values[index]

  # Add observer time
  def add_observer_time(self, observational_time:float, units:str):
    self.observer['observational_time'] = {}
    self.observer['observational_time']['unit'] = units
    self.observer['observational_time']['value'] = observational_time

  # Add system size
  def add_system_size(self, system_size:float, units:str):
    self.observer['system_size'] = {}
    self.observer['system_size']['unit'] = units
    self.observer['system_size']['value'] = system_size


################
## Function that generates the system of differential equations given a model
################
def generate_model_system(model:dict):
  """
  Function generating a function for the differential equations of a given model.

  Args:
      model (dict): Dictionary describing the model

  Returns:
      str: Python code defining the function for the model
  """
  system = ''

  # Iterate over the dependent variables in the model
  for variable in model['variables']:

    # Load dictionary associated with the variable
    variable_dict = model['variables'][variable]

    # Chaeck is independency and accept if False
    independent = variable_dict['independent']
    if not independent:

      # Build equation for the variable
      symbol = variable_dict['symbol']
      equation = symbol+'_dot = '
      for mechanism in variable_dict['mechanisms']:
        expression = variable_dict['mechanisms'][mechanism]
        equation += expression+' '

      # Update system
      system += equation+'\n'
  return system


########
## Function to obtain the symbols of the variables
########
def get_variables_symbols(model:dict):
  """
  Function that obtains the symbols of the variables of a model.
  Uses a model as input and returns a tuple with the symbols of the variables.

  Args:
      model (dict): Dictionary describing the model

  Returns:
      Tuple: Two objects
        Str: Symbol for the independent variable
        List: Symbols for the dependent variables
  """

  # Set list for dependent variables
  dependent_variables_symbol = []

  # Obtain variable symbols
  for variable in model['variables']:

    # Obtain dictionary for each variable
    variable_dict = model['variables'][variable]

    # Distinguish between independent and dependent variables
    independence = variable_dict['independent']
    if independence:
      # Obtain and store symbol
      variable_symbol = variable_dict['symbol']
      independent_variable_symbol = variable_symbol

    else:
      # Obtain and store symbol
      variable_symbol = variable_dict['symbol']

      # Append symbol and value to list
      dependent_variables_symbol.append(variable_symbol)

  # Return variables
  output = (independent_variable_symbol,dependent_variables_symbol)
  return output

################
## Function to obtain the parameter values of the model
################
def get_model_parameters(model:dict):
  """
  Function that obtains the peramaters of a given model.
  It expects a model as a dictionary and returns the symbols and values each as a tuple.

  Args:
      model (dict): Dictionary describing the model

  Returns:
      List: Two tuples
        tuple: symbols associated with the parameters.
        tuple: values associated with the parameters.
  """

  # Obtain dictionary for the parameters
  parameters = model['parameters']

  # Set the lists that later will become tuples
  param_symbols_list = []
  param_values_list = []
  for parameter in parameters:
    # Obtain the dictionary for the specific parameter
    parameter_dict = parameters[parameter]

    # Obtain and store symbol
    symbol = parameter_dict['symbol']
    param_symbols_list.append(symbol)

    # Obtain and store value
    value = parameter_dict['value']
    param_values_list.append(value)

  # Print what was obtained
  print(f'Parameter symbols: {param_symbols_list}')
  print(f'Parameter values: {param_values_list}')

  # Return the results as tuples
  param_symbols = tuple(param_symbols_list)
  param_values = tuple(param_values_list)
  output = [param_symbols,param_values]
  return output



############
## Function that returns the initial conditions for the variables in the model
############
def get_initial_conditions(model:dict):
  """
  Function that returns the initial conditions of the variables.

  Args:
      model (dict): Dictionary describing the model

  Returns:
      Tuple: Two lists
        List: Symbols for the initial conditions
        List: Values for the initial conditions
  """

  # Initiate lists
  initial_conditions_list = []
  initial_variables_list = []

  # Obtain dictionary for variables and loop through them
  variables = model['variables']
  for variable in variables:
    # Obtain dictionary for each variable
    variable_dict = variables[variable]

    # Obtain and store symbol
    symbol = variable_dict['symbol']
    initial_variables_list.append(symbol+'0')

    # Obtain and store value
    value = variable_dict['initial_value']
    initial_conditions_list.append(value)

  # Print output
  print(f'Initial symbol: {initial_variables_list}')
  print(f'Initial values: {initial_conditions_list}')

  # Prepare output as tuple for return
  output = (initial_variables_list,initial_conditions_list)
  return output


#######
## Function to generate the model equation as a python function .py
#######
def generate_model_fun(param_symbols:tuple,depend_variables_symbol:list,model_eqs:str):
  """
  Function that generates the python function for the model in a .py file.
  It takes the parameter symbols, dependent variable sybmols, and system of equations.
  It returns a string with the model and generates the python function as a file.py
  !!!! If the model was defined as a class with methods, this function would be more straightforward.

  Args:
  """
  ## Initialize output string
  output = ''

  ## Open file
  name = 'model_fun'
  file = open(f'{name}.py','w')

  ## Write definition line
  params = ','.join(param_symbols)
  definition = f'def {name}(t,y,{params}):' + '\n'
  file.write(definition)
  output += definition

  ## Write initial values
  index = 0
  print(depend_variables_symbol)
  for symbol in depend_variables_symbol:
    print(symbol)
    line = '\t' + f'{symbol} = y[{index}]'  + '\n'
    output += line
    file.write(line)
    index += 1

  ## Write system of equations

  line = '\n'
  output += line
  file.write(line)
  eqs_list = model_eqs.split('\n')
  for equation in eqs_list:
    eq_line = '\t' + equation  + '\n'
    output += eq_line
    file.write(eq_line)

  ## Write formatted output for y_dot
  y_dot = 'y_dot = ['
  for equation in eqs_list:
    # Obtain and recast the symbols for the derivatives to the generic y_dot
    variable_dot = equation.split('=')[0]
    y_dot += variable_dot +','

  y_dot = y_dot.rstrip(',')
  y_dot += ']'
  y_dot_line = '\t'+ y_dot  + '\n'

  output += y_dot_line
  file.write(y_dot_line)

  ## Write return output
  return_line = '\t'+ 'return y_dot'
  output += return_line
  file.write(return_line)

  ## Close file
  file.close()

  ## Return output
  return output

# MAIN

#------------------------------------------------------------------------------------------------
# Set model
# User input

## Choose if the model will be uploaded
upload = False ## Values: `True` or `False`. If True, it expects to upload model as JSON file.

## If upload is False, edit the lines below to define the new model
### Set variables info
variables=['Time','w_daysy','b_daysy']
variables_symbols = ['t','w_daysy','b_daysy']
variables_units = ['h','m^2','m^2']
variables_independence = [True,False,False]
variables_with_mechanisms = variables[1:]
mechanisms = [['growth_w','competition_w_b','carrying_capacity_w','temp_corr_growth_w','temp_corr_competition_w','temp_corr_carrying_w','decay_w'],
              ['growth_b','competition_b_w','carrying_capacity_b','temp_corr_growth_b','temp_corr_competition_b','temp_corr_carrying_b','decay_b']]
expressions = [\
['beta_max*p*w_daysy','-beta_max*b_daysy*w_daysy','-beta_max*w_daysy*w_daysy',\
'-p*k*(T-T_opt)**2*w_daysy', '-k*(T-T_opt)**2*b_daysy*w_daysy', '-k*(T-T_opt)**2*w_daysy*w_daysy',\
'-gamma*w_daysy'],
['beta_max*p*b_daysy','-beta_max*w_daysy*b_daysy','-beta_max*b_daysy*b_daysy',\
'-p*k*(T-T_opt)**2*b_daysy', '-k*(T-T_opt)**2*w_daysy*b_daysy','-k*(T-T_opt)**2*b_daysy*b_daysy',\
'-gamma*b_daysy']]
variables_values = [0,0.1,0.2]


## If upload is False, edit the lines below to define the new model
### Set variables info
# variables = ['Time','Bacteria','Phage']
# variables_symbols = ['t','B','P']
# variables_units = ['h','cells/ml','phages/ml']
# variables_independence = [True,False,False]
# variables_with_mechanisms = variables[1:]
# mechanisms = [['growth','predation'],
#                 ['burst','decay']]
# expressions = [['+r*B','-aB*B*P'],
#                  ['c*aP*B*P','-m*P']]
# variables_values = [0,100,100]

### Set parameters info
parameters = ['growth_rate','free_soil','temp_constraint',
                'temp','optimal_temp','decay_constant']
parameters_symbols = ['beta_max','p','k','T','T_opt','gamma']
parameters_units = ['1/day','m^2', '1/K**2','K','K', '1/day']
parameters_values = [1, 1, 0, 280, 295.5, 0.05]
#parameters_values = [1, 1, 17.5**(-2), 280, 295.5, 0.5]


### Set parameters info
# parameters = ['growth_constant','infection_constant_to_B',
#                 'burst_size','infection_constant_from_P','decay_constant']
# parameters_symbols = ['r','aB',
#                         'c','aP','m']
# parameters_units = ['1/h','ml/(h*phage)',
#                       '','ml/(h*cell)','1/h']
# parameters_values = [1,2,3,1,2]

### Set observer info
observational_time = 100
observational_time_units = variables_units[0]
system_size = 1
system_size_units = 'm^2'
              
### Set observer info
# observational_time = 10
# observational_time_units = variables_units[0]
# system_size = 1
# system_size_units = 'ml'


## If upload is True the coad below will request to upload the model (no need to modify)
if upload:
  # Upload model from JSON file
  uploaded = files.upload(model_fun)
  model_bytes = next(iter(uploaded.values()))
  model = json.loads(model_bytes.decode())
  print('Model uploaded')


## The code below will run to generate the model if the user chose upload=False.
if not upload:
  # Set model instance
  model_object = Model()

  ## Set observer info
  observational_time = 100
  observational_time_units = variables_units[0]
  system_size = 1
  system_size_units = 'ml'

  ## Load variables info to the model
  model_object.add_variables(variables)
  model_object.add_symbols_to_variables(variables,variables_symbols)
  model_object.add_units_to_variables(variables,variables_units)
  model_object.add_status_to_variables(variables,variables_independence)
  model_object.add_mechanisms(variables_with_mechanisms,mechanisms,expressions)
  model_object.add_initial_values_to_variables(variables,variables_values)

  ## Load parameters info to the model
  model_object.add_parameters(parameters)
  model_object.add_symbols_to_parameters(parameters,parameters_symbols)
  model_object.add_units_to_parameters(parameters,parameters_units)
  model_object.add_initial_values_to_parameters(parameters,parameters_values)

  ## Load observer info to the model
  model_object.add_observer_time(observational_time,observational_time_units)
  model_object.add_system_size(system_size,system_size_units)

  model = model_object.__dict__

  ## Store model
  with open(f'model.json', 'w') as file:
    json.dump(model, file, indent=4)
#------------------------------------------------------------------------------------------------

# Display model
print(model)

# Generate and display equations
model_eqs = generate_model_system(model=model)
print(model_eqs)

# Load the model as a function that can be integrated

## Get symbols for dependent variables and
indep_variable_symbol, depend_variables_symbol = get_variables_symbols(model)
param_symbols, param_values = get_model_parameters(model)
model_fun_str = generate_model_fun(param_symbols,depend_variables_symbol,model_eqs)
print(model_fun_str)
exec(model_fun_str)

# Conditions

## Parameters
[param_symbols,param_values] = get_model_parameters(model=model)

## Initial values
inital_symbols, initial_values = get_initial_conditions(model=model)
y0 = initial_values[1:]

## Time conditions
## Obtain time conditions from model
t0 = model['variables']['Time']['initial_value']
tend = model['observer']['observational_time']['value']
tspan = [t0,tend]
max_step = 0.01
print(f't0 = {t0}')
print(f'tend = {tend}')
print(f'max_step = {max_step}')

args = param_values
#solution = solve_ivp(test_model, tspan, y0, args=args, dense_output=True,max_step=max_step)
fun = model_fun
solution = solve_ivp(fun, tspan, y0, args=args, dense_output=True,max_step=max_step)

# Plot phase diagram dynamics
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
title = 'Phase diagram'
x = solution.y[0]
y = solution.y[1]
label_x = depend_variables_symbol[0]
label_y = depend_variables_symbol[1]
plt.plot(x, y, 'b-')
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel(label_x)
plt.ylabel(label_y)
plt.title(title)
plt.legend()
plt.show()

# Plot daysies dynamics
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
title = 'Daysies'

time=np.linspace(t0,tend,len(x))

x = solution.y[0]
y = solution.y[1]

label_x = 'time (day)';label_y = 'daysy cover (m$^2$)'
plt.plot(time, x  ,'limegreen'   , label='clear')
plt.plot(time, y  ,'navy', label='dark')
plt.plot(time, y+x,'k-'   , label='total')
y_pos=[0, 0.25, 0.5, 0.75, 1]
ax.set_yticks(y_pos)
plt.xlabel(label_x)
plt.ylabel(label_y)
plt.title(title)
plt.legend()
plt.show()


