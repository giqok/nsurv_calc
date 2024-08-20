# Calculating monotonically decreasing survival from data with missing observations

## Set up for the data
The data needs to be structured with a numeric 'time' column and columns for each population compartment as its own column. 
Each population compartment must be mututally exclusive (i.e. the sum across all compartments must be less or equal to the total initial population).
Treatments and Replicates need to be noted as such, either as categorical columns or as numeric concentrations.
Note that replicates mean within treatment replicates, and not an integer representation of treatment number.

## Set up for the function
Provide the column names for observed survivors, these will be summed together to a grand total number of survivors in 'survivor_columns = '. 
Likewise, provide the column names for observed dead (or immobilized in case of non-recovery) in 'dead_columns = '.
In case of a single survival or dead column, provide a simple string, otherwise a list of strings.
If there is a compartment sink for individuals that have grown out of or entered into that you'd like to not be counted, provide a list of these column names in 'sink_columns = '.
Set 'initial_population = ' to the column name of initial population counts in case you want to assume the first count is identical to the number of individuals introduced at the beginning of the experiment, otherwise this is estimated from the maximum known living count (see assumption one).

## Function output
An additional column for Nsurv is provided for each treatment:replicate:time combination.
This column is monotonically decreasing across time for each treatment:replicate.

## Function details
First assumption, no zombies.
If a number of surviving members is found at a later time, no fewer than that amount can be alive at an earlier time point. 
The known living count is the max of living counts for all time points from now to the last observation time. 

Second assumption, the difference between known living and initial population defines the max total number of potential dead.
This may contradict the 'dead' count observed at a given time point, but the first assumption implies that the dead count may include a reversed immobilization. 

Third assumption, counts of dead individuals may not be exclusive of individuals already counted.
This is generally true when the dead individuals are not removed from the system.
The second assumption defines an upper limit on how many dead we accept in calculating the final Nsurv count, and the cumulative sum of dead over the current and all previous observations is an upper limit on the observed dead individuals.
This method also covers the case when the dead are removed or otherwise not counted in future observations.
The expected dead count is the minimum of cumulative observed dead and maximum total number of potential dead as defined in the second assumption.

As a corollary, the number of known living and number of expected dead together is never higher than the initial population.

The final Nsurv is calculated as initial population minus the expected dead count at every time point. 

## The math
For any given treatment:replicate, there are time points $t_i, i=0,1,\dots,n$ for every observation.
Survivors observed at time $t_i$, $S_i=S(t_i)$, is the sum of all compartments observed across the 'survivor_columns'.
Initial survivors, $S_0=S(t_0)$ is either given by the column 'initial_population' or taken as the max $S_0=\max_{i}SK(t_i)$, where $SK$ is defined below.
Similarly, the dead observed, $D_i$, is the sum of all compartments observed across the 'dead_columns'.
The known survivors, $SK$, at time $t_j$ is calculated as

$$SK(t_j)=max_{j\leq i \leq n}S_i$$

and the maximum potential dead, $DP$ at time $t_j$ is calculated as

$$DP(t_j)=S_0-SK(t_j).$$

This becomes the upper bound we compare the cumulative sum of observed dead, calculated as

$$DOS(t_j) = \sum_{i=0}^{j}D_i$$

which gives us the expected dead

$$DE(t_j) = \min(DOS(t_j), DP(t_j)).$$

The final survival count, $Nsurv$, at each time $t_j$ is calculated by 

$$Nsurv(t_j) = S_0-DE(t_j)$$

which is then returned as a new column on the given data.frame.

## The code
```
library(dplyr)

nsurv_calc <- function(df,
                       treatment,
                       replicate,
                       time,
                       survivor_columns,
                       dead_columns,
                       initial_population) {
  out_df <- df |>  # Rowwise summation across provided columns
    mutate(total_surv = rowSums(df |>
                                  select(all_of(survivor_columns)), na.rm = TRUE)) |>
    mutate(total_dead = rowSums(df |>
                                  select(all_of(dead_columns)), na.rm = TRUE)) |>
    arrange(treatment, replicate, time)  # order for cumulative calculations
  
  out_df <- out_df |>  # group and calculate cumulative values
    group_by(.data[[treatment]], .data[[replicate]]) |>
    mutate(known_survivors = rev(cummax(rev(total_surv)))) |>
    mutate(sum_of_dead = cumsum(total_dead))
  
  if (exists('initial_population')) {
    # use an initial population column
    out_df <- out_df |>
      mutate(max_possible_dead =
               pmax(0, .data[[initial_population]] - known_survivors)) |>
      mutate(expected_dead = pmin(max_possible_dead, sum_of_dead)) |>
      mutate(Nsurv = .data[[initial_population]] - expected_dead)
  } else {
    # Otherwise use max known_survivors
    out_df <- out_df |>
      mutate(max_possible_dead = max(known_survivors) - known_survivors) |>
      mutate(expected_dead = pmin(max_possible_dead, sum_of_dead)) |>
      mutate(Nsurv = max(known_survivors) - expected_dead)
  }
  
  return(
    out_df |>
      select(
        -expected_dead,
        -max_possible_dead,
        -sum_of_dead,
        -known_survivors,
        -total_surv,
        -total_dead
      )
  )
}
```
