# Single-cell fluorescence trace fitting
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1418465.svg)](https://doi.org/10.5281/zenodo.1418465)

This repository provides code for fitting model functions to single-cell
fluorescence time traces using multistart maximum-likelihood fitting.


## Installation
This program requires Matlab in a recent version. It was developed and
employed in R2016b and R2017a. However, other versions may work as well.

To install the program, clone the repository to a place accessible to
Matlab.

To ensure that all necessary scripts are on your Matlab path,
run the `install_apoptosis_fit` function.

Note that the files in this repository use the UTF-8 encoding.
If your copy of Matlab is not configured to use UTF-8, some characters
may be displayed incorrectly.
However, the program should still operate as intended in most cases.


## Usage
The program can be run interactively or as a function.

The interactive mode is useful for a small number of data files and can be
used without much knowledge. However, the interactive mode only allows
using one model per run and saves the output files to a default directory.

Running the program as a function provides more flexible options.
For example, it is possible to process a large number of files spread over
multiple directores, employ multiple models in one run and to specify a
custom directory for result output.


### Data prerequisites
The files containing data to be fitted must be CSV or XLS(X) files.

In both cases, only numeric data and no header lines are allowed.
The first column must contain the time values of the data points.
The other columns are the traces in the file, so that the (n+1)-th column is the n-th trace in the file.
The time should be given in hours. Refer to the description of the `TimeUnit`, `TimeScale` and `TimeOffset` options for details.
The fluorescence values can be given in arbitrary units.

CSV is the preferred format.
The input CSV files must use the default Matlab formatting, which is a `.` as decimal separator, a `,` as field/column separator and a newline character as line separator.

XLS(X) has not been used and tested as thoroughly as CSV.
Moreover, the behaviour may be platform specific.

There is [example data](https://doi.org/10.5281/zenodo.1418377) available.
To use it, download one of the ZIP files and extract it.
The example data for fitting is below the directory `Raw`.
Note, however, that when using example data from `Data_Huh7.zip`, the value of `min_breakdown_amp` in [postproc_kink.m](postprocessing/postproc_kink.m) should be changed in order to account for the different behaviour of Huh7 cells:
```diff
- min_breakdown_amp = .1 * R.amplitude;
+ min_breakdown_amp = .2 * R.amplitude;
```


### Interactive mode
The interactive mode requires that Matlab can display figure windows.

Simply execute [fit_apoptosis](fit_apoptosis.m), for example by typing
`fit_apoptosis` in the Matlab command window.

A dialog window appears that asks you to select a model. Select the model
you want to use for fitting your data and click “OK”.

Then, a dialog window appears that asks you to select the data files
containing the traces you want to fit. To select multiple files, you can
press Ctrl while clicking on the file. The selected files must contain the
data as described in “Data prerequisites”. After selecting the files,
click “Open”.

Then the files are read in and processed automatically.

If you have Matlab’s “Parallel Computing Toolbox” and a corresponding
license installed, a parallel pool is created using the default profile,
if not running yet, and the fitting is performed in parallel.
To use a custom parallel pool profile, start the parallel pool before
invoking `fit_apoptosis`.
If no parallel pool can be started, the fitting is performed serially.

When the fitting is finished, the files described in “Output files” are
created in the directory `project/results/<model>/<date>`, where
`<model>` is the name of the model you selected, and `<date>` is the
date at which you have started the fitting.


### Function call
Calling `fit_apoptosis` as a function is especially suitable for usage
in a script, where no interaction with the user is desired.

The following function signatures are supported:
```
fit_apoptosis()
fit_apoptosis(model)
fit_apoptosis(model, file_paths)
fit_apoptosis(model, file_paths, out_path)
fit_apoptosis(model, file_paths, out_path, run_options)
fit_apoptosis(model, file_paths, out_path, run_options, export_options)
```
* `model` is the model to be fitted. It is either a char array (one model)
   or a cellstring of one or more model names. If only one model is
   specified, it is used for all data files. If multiple models are
   specified, the number of models must match the number of data files, and
   each data file is fitted with the corresponding model. If no model is
   specified, a model selection dialog is shown.

* `file_paths` are the paths of the files containing the traces to be
   fitted, as described in “Data prerequisites“. It is either a char array
   (single data file) or a cellstring. If no file is given, a file
   selection dialog is shown.

* `out_path` is the directory to which the output files are written. It is
   specified as a char array (single output directory) or a cellstring.
   If multiple output directories are specified, their number must match the
   number of datafiles given.
   If multiple models are specified, there must also be specified multiple
   output directories.
   If no output directory is specified, the default output directory as
   described in “Interactive mode” is used.

* `run_options` are options for fine-tuning the behaviour. It is a cell
   array of keys and values (`{'key1', value1, 'key2', value2}`). The keys
   are char arrays. The following keys are supported:
   * `TempDir`: the path of a directory to be used for saving the temporary
      files, given as cellstring. Default is the output directory.
   * `TimeUnit`: the char array to be used for the time unit in the plots.
      Default is “h” for hours.
   * `TimeScale`: the value with which to multiply the time values read
      from the data files to convert them into hours. By default, the time
	  scaling is recognized by the function `get_t_scale` in the file
	  [private/scale_t.m](private/scale_t.m).
   * `TimeOffset`: a value added to the time values, in hours. Default is 0.

* `export_options` are options for controlling what output files are
   exported. By default, all possible output files are generated.
   If `export_options` is a char array, of a option name, only the
   corresponding option is enabled; all others are disabled. The value
   `'all'` enables all options. If `export_options` is a cellstring, the
   options included in it are enabled and all others disabled. If
   `export_options` is a structure with option names as fields, the default
   values are overwritten with the values in the structure.

   The possible options are:
   * `'debug'`: plot detailed single-trace graphs with model-specific
      lines and visualizations.
   * `'single'`: plot single-trace graphs showing raw and fitted traces.
   * `'total'`: plot an overview graph for each data file, containing all
      raw traces in the data file and the corresponding fits.
   * `'params'`: export a CSV table containing the parameter best
      parameter estimates for all traces.
   * `'states'`: export a CSV table containing postprocessing and filtering
      results for all traces.
   * `'simulated'`: export CSV files containing fits of the traces.
   * `'temp'`: leave the main temporary file after finishing, instead of
      deleting it.

   The structure of the content of the output files is described in
   “Output files”.


### Batch processing
To leverage batch processing, the fitting can be performed separately
from exporting the data. For this purpose, the function
`fit_apoptosis_batch` can be used. It takes the same arguments as
`fit_apoptosis` and stops after fitting and postprocessing is finished.

`fit_apoptosis_batch` creates the directory hierarchy for the output files
and a file called `TEMP_<hhmm>.mat`, where `<hhmm>` is the time in hours
and minutes when the fitting was started. The timestamp is added in order
to avoid overwriting temporary files from other runs in the same directory.
Fitting and (parts of) postprocessing are performed in a `parfor` loop
under the same conditions as described for the interactive mode.

To create the output files, call `fit_apoptosis_output` with the path to
the `TEMP_<hhmm>.mat` file as argument. Then, the prepared output directory
will be populated with the output files.

Note that the output
directories are saved as absolute paths. When copying the temporary file
to another computer for plotting, the path names must be the same on both
computers. For Linux users, it is recommended to use the tilde to denote the home
directory of the current user for specifying the output directory.

Since sporadic Matlab crashes were encountered during plotting, the index
of the trace or file whose data are currently exported is written to the temporary file `plot.log`.
From there, it is automatically recovered when `fit_apoptosis_output` is called
again after a crash, and plotting/exporting is resumed where it had stopped before.

However, if no `plot.log` file is found but the output directories contain CSV or PDF files,
an error is thrown since you might be overwriting already exported files.
In this case, the index of the trace or file at which exporting should be
started must be specified manually by setting `loop_start`,
the second (and optional) parameter of `fit_apoptosis_output`,
to the trace number indicated in `plot.log`.
If `plot.log` indicates a file number, `loop_start` should be a vector with `NaN`
as its first element to skip exporting single trace data, and with the
trace number as its second element.
(The single trace data are exported before file data.
Hence, if Matlab crashes during exporting file data, the single trace data
has already been exported an needn’t be re-exported.)


### Output files
`fit_apoptosis` and `fit_apoptosis_output` write some files to the output
directories. Which files are exported, can be controlled by specifying
the `export_options` argument, as described above. By default, all types
of output files are exported.

The output files have names according to a scheme:

    <name>_<TYPE>_[<number>_]<hhmm>.<suffix>

* `<name>` is the file name of the corresponding input file without suffix.
* `<TYPE>` indiates which kind of data is exported to this file and is written in capital letters.
* `<number>` is only present in the names of output files corresponding to a single trace instead of all traces in the data file. It is the number of the trace in the file, starting at 1, possibly zero-padded. (The square brackets are not part of the file name; they only indicate the part that is not always present.)
* `<hhmm>` is a timestamp of the time when the fitting was started. It is used to prevent overwriting files when the same files are fitted again.
* `<suffix>` is the file suffix, which indicates the file type. For data files that are meant for further processing, it is `csv`. These are valid CSV files in Matlab’s default format: The decimal separator is `.`, the field/column separator is `,`, and the line separator is a newline character. These files only contain numeric data. For files containing plots, the suffix is `pdf`. These are PDF files generated by Matlab.

The following values for `<TYPE>` are possible:

`ALL_FIT` is a PDF file contining an overview plot.
It shows all traces in a file and their corresponding fits.
The measured data are shown as thin lines, the fits as thick lines.
Different traces are shown in different colors.
Since, however, only a limited set of colors is available, the same color may be used for several traces.
The overview plot may be used to compare the traces to each other, to find outliers or to get a general impression of the traces in the corresponding file.

`ALL_PARAMS` is a CSV file containing the estimated values of the parameters for each trace.
The line index corresponds to the trace index in the file, and the column index corresponds to the parameter index, where the parameter order defined for the model is used.

`ALL_SIMULATED` is a CSV file containing fitted traces.
It may be used for plotting the fits without re-evaluating the model function.
The first column contains the time values in hours.
All other columns correspond to traces, so that the (n+1)-th column corresponds to the n-th trace.
The lines correspond to different time points.
The time resolution is currently ten times the time resolution of the original data file.

`ALL_STATE` is a CSV file containing information obtained by the postprocessing procedure.
* The first column is the index of the trace in the file (one-based).
* The second column is the event time recognized; failure of event time recognition is indicated by conventional values, typically `NaN`, `Inf` or `-Inf`, which may code a reason for the failure.
* The third column contains the absolute amplitude of the trace, i.e. the difference between the largest and the smallest data point in the trace.
* The fourth column contains the relative amplitude of the trace, i.e. the absolute amplitude divided by the amplitude of the trace with the largest non-outlier amplitude in the file.
* The fifth column contains the logarithmic likelihood of the best fit.
* The sixth column contains the fit type.
   This is a numeric identifier of the fit function used and of the postprocessing routine used.
   The following values are currently defined:
   * `1`: Indicates postprocessing by `postproc_LATE_extrapol.m`.
   * `2`/`-2`: Indicates postprocessing by `postproc_kink.m` or `postproc_tmrm_kink.m`.
   * `3`: Indicates postprocessing by `postproc_mRNA_trivial.m`.
   * `4`: Indicates postprocessing by `postproc_LATE_decay.m`.
   Custom models should define another, unique identifier to facilitate identifying the model used from within the data.
* The seventh column contains the slope of the trace at the event.

`DEBUG` is a PDF file containing a detailed plot of one trace.
It may contain marks defined by the postprocessing to visualize trace properties or the postprocessing results in a detailed way.
As the name suggests, this file can be used to analyze and improve the postprocessing algorithms.

`FIT` is a PDF file containing a plot of one trace.
The measured trace is drawn in blue and the fitted trace in red.
This file may be used for publishing single traces, for example.


## Implementing custom models
To implement custom models, the following changes are necessary:

1. Implement a model function to be fitted.
    Skip this step if you want to use an already existing model function.
    Model functions are conventionally stored in the `simulations` directory.
    The model function takes two arguments: `t`, which is a vector of times (in hours) at which the model function is to be evaluated, and `params`, which is a vector of values for the model parameters.
    It must return a column vector of as many elements as in the time vector, and each element must contain the value obtained from the model for the given parameters at the given time.

2. Implement a postprocessing function.
    Skip this step if you want to use an already existing postprocessing function.
    The postprocessing function takes three arguments: `D`, `F` and `R`.
    Each of the arguments is a structure.

    `D` contains data about the raw trace. The field `data` is a vector of the raw trace (to which a additive offset may be applied to eliminate negative values in the corresponding file), and the field `index` is the index of the trace in the file.

    `F` contains information about the file in which the trace is contained.
    The field `t_sim` contains a vector of times with a ten-fold higher resolution than in the original data file.
    The field `modelname` the name of the model used for processing the file as an char vector, and the field `simulate` contains a file handle to the model function.

    `R` contains fitting results.
    The field `params` is the vector of the best estimates for the model parameters.
    The field `data_sim` is a vector containing the values obtained by evaluating the model function with the parameters from `R.params` at the times in `F.t_sim` and may be used as a highly resolved representation of the fitted trace.
    The field `min_val` is the smallest value in `data_sim`.
    The field `min_ind` is the index of the smallest value in `data_sim`.
    The field `max_val` is the largest value in `data_sim`.
    The field `max_ind` is the index of the largest value in `data_sim`.
    The field `amplitude` is the amplitude of the fit, which is the difference between `max_val` and `min_val`.
    The field `data_amp` is the difference between the largest and the smallest value in `D.data`.

    The postprocessing function should return a `varargout`.
    The first output should be the event time recognized.
    The second output should be a struct comprising model dependent information. The following fields are available, which may have a scalar value: `parabola`, `rising_edge_xoff`, `rising_edge_yoff`, `rising_edge_scale`, `falling_edge_xoff`, `falling_edge_yoff`, and `falling_edge_scale`.
    Moreover, the field `type` should contain a numeric identifier of the model used; it is the same that is written to the sixth column of the STATE output file.
    The third output can be a model dependent scalar value such as the slope of the trace at the detected event.

    The postprocessing files are, by convention, located in the directory `postprocessing`.
    You may want to have a look at the existing postprocessing files to get some inspiration on how a postprocessing file may be structured.

3. Set up the debug plotting.
    If you want to have a verbous DEBUG plot file, you have to add instructions specific for your model to the function [plot_postproc](postprocessing/plot_postproc.m).
    You can orientate yourself by the existing content of this function.

4. Define your actual model.
    Create a new entry for your model in the function [model_definitions](model_definitions.m).
    Your model simply needs its own `elseif` condition block.
    The condition must be true if the variable `modelname` equals the name of your model.
    This can be done by using `strcmp`.
    However, to make your model more robust against typos or to allow for abbreviations of the model name, you may want to perform a case-insensitive comparison or even use regular expressions.

    Then, you only have to fill in some model-specific data into the fields of the struct `model`, just as you can see for the existing models.
    The field `name` is the name of your model that will be displayed for identifying your model to the program user.
    The field `marker` is the name of the marker used in your measurements; it will be written into the plots. The marker may or may not be equal to the model name.
    The field `simulate` is a function handle of the model function you implemented (or chose) before.
    The field `postproc` is a function handle of the postprocessing function you implemented (or chose) before.
    The field `par_num` is the number of model parameters used for fitting.
    Note that the last parameter will be used as variance of the residual distribution for fitting, so you have to provide one more parameter than your model requires.
    The fields `par_min` and `par_max` are vectors defining the logarithms of the smallest and the largest value, respectively, the parameters can take.
    The lengths of these vectors must be equal to the number of parameters.
    Note that the values in the vectors are (decadic) logarithms of the values at which your model function is evaluated.
    The field `par_names` is a cell vector of parameter names which may eventually be used for detailed plotting.
    You can denote mathematical expressions in simple TeX-like syntax.

5. To allow for selection of your custom model in an interactive session, add your model name (that is matched by the `elseif` condition) to the list of models in line 94 in the function [readInputData](private/readInputData.m).

That’s all.
Now your custom model is readily set up for use.

