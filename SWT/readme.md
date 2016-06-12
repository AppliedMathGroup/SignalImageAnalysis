
MATLAB function using the Square Wave Method (SWM) for time series and signal analysis
==================

####Ricardo E. Monge ,   Sherry Gapper,    Osvaldo Skliar

The MATLAB® command described here performs a time series analysis using the Square Wave Transform (SWT)<sup>1</sup>, and is accessible from the MATLAB interface under the name `swt`.

`swt` is a user-defined function that requires three parameters to work properly. The first parameter `V` corresponds to the time series to be analyzed, input as a standard MATLAB vector. The second parameter `f` corresponds to the sampling frequency in Hz of the data in parameter `V`. Finally, the third parameter `Dt` corresponds to the time interval in seconds of the entire data set. 

The outputs are a two-column matrix `S`, and a single numerical value that indicates the approximation quality, referred to as `Dm`. The matrix `S` will have a row for each data value in the input data set, in which the value in the first column correspond to the square wave frequency (f<sub>i</sub> in the published paper) and the value in the second column corresponds to the square wave coefficient (C<sub>i</sub> in the published paper).

If, within a certain time interval, a given time series is an adequate representation of a signal, then the SWM analysis of that time series, in that same time interval, can also be considered an analysis of that signal in that interval.

The approximation obtained with the SWM for each of the numerical values in the time series analyzed is outstanding. Let V<sub>i</sub> (for i = 1,2,…,n) be the *i* th value of the time series analyzed, and V<sub>icomp</sub>, the corresponding approximation computed by adding all the square waves (whose coefficients and frequencies were computed by the SWM implementation) for that point. Thus, |V<sub>i</sub> - V<sub>icomp</sub>| is the modulus (absolute value) of the difference between those two values. There will be a modulus of this type for each data point in the time series. We can, therefore, compute a maximum of this difference between the time series:

D<sub>m</sub>=max |V<sub>i</sub> - V<sub>icomp</sub>|.

Therefore, each computation using the `swt` function defined here will have the following outputs:

1. The sequence of dyads, in the form of “frequencies - coefficients” (f<sub>i</sub> - C<sub>i</sub>) for each of the values in the time series being analyzed. This sequence provides a detailed view of the time series in the frequency domain. This representation is known as Square Wave Transform (SWT) and its acronym is used for the MATLAB command (`swt`). 
2. The computed value for D<sub>m</sub>.

####References
1. [Skliar O., Monge R. E., Oviedo G. and Gapper S. (2016). A New Method for the Analysis of Signals: The Square Wave Transform, *Revista de Matemática. Teoría y aplicaciones* 2016, Vol. 23(1), pp. 85-110.](http://revistas.ucr.ac.cr/index.php/matematica/article/download/22352/22509)

