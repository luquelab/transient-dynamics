# WORKFLOW FOR HYPERION

## Step 1 - Characterize the system
The first step in studying an ecosystem (or any other physical system, for that matter) is to know what its constituting elements are: plants, animals, bacteria, phages...For Hyperion, these are essentially the **Things**.  For instance, in the consumer-resource system we are considering, we have two **Things**: phages and
sensitive bacteria. How do **Things** interact? **Things** interact through **processes**.  For instance, in our system we have four **processess**: growth,
infection, burst and viral decay.


## Step 2 - Modeling
Considering the **Things** and **Processess* found on Step 1, we build a mechanistic dynamic model. That is, a model that considers time explicitly and explicitly
includes the mechanisms that control the system. In the system we are considering, Lotka-Volterra equations model the behavior of the system.

## Step 3 - Critical concentrations
Mechanistic models have different terms that compete and interact with each other. Likewise, the contribution of these terms to the dynamics of the system changes
with time and depending on the specific parameters. More specifically, these changes occur when critical concentrations are reached.
Critical concentrations can be defined as the concentrations that effectively inactivate per-capita **processes**. Per-capita **processes** are simply how processes impact single bacterial or phage entities. Take again our model, for example. If the concentration of bacteria is of the order of 'K', the carrying capacity, that inactivates the growth term of the model. On the other hand, if the phage concentration is close to r/d (r and d being the growth and infection rates, respectively) the infection part of the model becomes inactive. Hyperion can actually track and control when critical concentrations are reached.

## Step 4 - History Dynamics
Critical concentrations can be estimated *in silico*. Three things are required:
1. A dynamic model
2. Numerical values for the parameters of the model. Alternatively, a reasonable range of values can be explored using Latin Hypercube Sampling, for instance.
3. Hyperion: Hyperion tracks the per-capita contribution of each **process** to the dynamics of the model-

The combination of these three elements outputs concentrations or abundances of **things** plus the time points at which they are reached. These time points are the target of statics.

## Step 5 - Statics Prediction
Using the outputs of the previous step, we can estimate the time snapshots (statics) at which experimental samples should be collected to gain the most insights into them. It is important to note that experiments can only measure **things** and not **processes**. Hyperion helps us to spot the needle in the haystack: what are the concentrations of **things** that reveals most information about the **processess**.