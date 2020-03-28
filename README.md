# Shockley Queisser limit : 
## Theoretical Dye-sensitized Solar Cell Efficiency calculator and visualizer

It calculates the theoretical dye-sensitized solar cell parameters with options to change temperature, light intensity, and radiative efficiency, provides visualization tools in python and MS Excel. 

The calculation is based on the radiative limit or the Shockley Queisser limit as described here [Wikipedia: Shockley–Queisser limit] (https://en.wikipedia.org/wiki/Shockley–Queisser_limit); however, we modified the calculation to estimate the theoretical efficiency limit of dye-sensitized solar cells.

This python script and the MS Excel sheet are based on the work by [C. Marcus Chuang](https://github.com/picklu/Shockley-Queisser-limit)



#  Example Outputs
## Bandap vs Efficiencies of a single junction solar cell

Create a `SQlim` object, which calculates all the parameters and then call the method `plot('PCE')` to generate the bandgap vs efficiency plot. PCE: power conversion efficiency.
```python
SQ = SQlim()
SQ.plot('PCE')
```

<img src="/ExampleOutputFig/PCE.png" width="600">


The efficiencies will be in the `SQ.PCE` attribute, an numpy array. There's also a dictionary `SQ.paras` that stores all the characteristics as a key (string): characteristics (numpy array) pairs.


## Four important parameters in a single figure: 

Calling the method `plotall()` will generate a plot containing the 4 most important characteristics of a solar cell in the subplots
```python
SQ.plotall()
```
* Voc: Open-circuit voltage,
* Jsc: Short-circuit current density,
* FF: Fill factor

<img src="/ExampleOutputFig/ALL.png" width="800">

A method `get_paras(self, Eg, toPrint=True)` can be used to look up the results. For example, the following call would print the theoretical parameters for a 1.337 eV solar cell.
```python
SQ.get_paras(Eg=1.337)
```
would print the following lines like these in the colsole:
```python
"""
Bandgap: 1.337 eV ; J0 = 2.64e-17 mA/cm^2

Voc = 1.079      V
Jsc = 35.14      mA/cm^2
FF  = 88.88      %
PCE = 33.703     %
"""
```


## Plot other characteristics

The `plot(para)` method can be used to generate different plots. Valid input `para` are `"Voc"`, `"Jsc"`, ``"FF"`, `"PCE"`, and `J0` (dark saturation current)
```python
SQ.plot('J0') # dark saturation current
SQ.plot('Voc') 
```

<img src="/ExampleOutputFig/J0.png" width="450"> 
<img src="/ExampleOutputFig/Voc.png" width="450">


# Calculate and plot the J-V curves

The `simulate_JV()` method can be used to calculate the J-V curve with an option to plot it.

```python
SQ.simulate_JV(Eg=1.337, plot=True) # calculate the J-V curve of solar cell, plot it, and return the J-V data
```

<img src="/ExampleOutputFig/JVcurve_1pt337eV.png" width="450">


##  Savedata
The data (Voc, Jsc, FF, PCE, J0 as a function of bandgap) can be saved as a single .csv file
```python
SQ.saveall(savename="SQ lim") # save data as "SQ lim.csv"
```

##  The data can be accessed here: [SQ limit data](/SQ%20limit.csv)

#
#
# Visualize more interesting results

## Break down of the theoretical efficiency and the energy loss

The method `SQ.E_loss(Eg)`, which takes bandgap `Eg` as an input, can be used to visualize the break down of energy loss. 

```python
SQ.E_loss(Eg=1.337)
```

Shown here are the break down for a 1.337 eV solar cell, which has the maximum theoretical efficiency of 33.7 %.

<img src="/ExampleOutputFig/E_loss_1pt337eV.png" width="800">


##  Available Energies

The mathod SQ.available_E(Egs) can be used to calculate and plot theoretical maximum available energies from a series of (mechanically stacked) solar cells with different Egs.

### Single-junction solar cell, 1.337 eV

```python
SQ.available_E(Egs=1.337)
```

This is the similar to the one above but without the break down of energy loss.

<img src="/ExampleOutputFig/E_avail_1pt337eV.png" width="800">

  
#
#
#
# Different Conditions

The default conditions for calculating the theoretical limits are the standard conditions : Temperature `T = 300` K, 1-sun condition `intensity = 1.0`, and radiative efficiency `EQE_EL = 1.0` (100% external quantum efficiency for electroluminescence). 
```python
class SQlim(object):
    def __init__(self, T=300, EQE_EL=1.0, intensity=1.0):
        """
        T: temperature in K
        EQE_EL: radiative efficiency (EL quantum yield)
        intensity: light concentration, 1.0 = one Sun, 100 mW/cm^2
        """
```

We can calculate the efficiencies in different conditions by simply changing the input. 

Because of this flexibility, we can easily get an idea of how the change of these factors affect the theoretical efficiencies (or other characteristics fo interest).

## Different Temperature

#### The *theoretical* maximum possible efficiencies of solar cells could be higher at lower temperature.
The classmethod `vary_temp()` can do this calculation and plot the results.

```python
SQlim.vary_temp(T=[150, 200, 250, 300, 350, 400])
```

<img src="/ExampleOutputFig/VaryT.png" width="800">


## Different Light intensity (Solar concentrator)

####  The efficiencies are higher when the incident light intensity is higher. 
The classmethod `vary_suns` does that caculation. 
```python
SQlim.vary_suns(Suns=[1, 10, 100, 1000])
```

<img src="/ExampleOutputFig/VaryIntensity.png" width="800">


## Different radiative efficiency

####  The lower the EQE_EL (radiative efficiency), the lower the power conversion efficiency.
```python
SQlim.vary_EQE_EL(EQE_EL=[1, 1E-2, 1E-4, 1E-6])
```
<img src="/ExampleOutputFig/VaryEQEEL.png" width="800">







