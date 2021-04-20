# Analysis
![status](https://github.com/TING2938/Analysis/workflows/C/C++%20CI/badge.svg)

**Provide some useful class and function to process the data from MD.** 

---

## class and function

- container of data
	- **Vector**
	- **Matrix**
	- **Tensor**
	
- **Getopt** : tools for parsing command line options.

- **GmxHandle** : trajectory analysis handle for the **Gromacs&reg;**  

---

- **loadtxt** : load ascii data to **Matrix**.
- **localTime** : get date and time.
- **arange, linspace, sum, max, min, mean, stdev, argMax, argMin, trapz, cumTrapz, cumSum, cumProduct, diff, print, find, swap, part, append, fill, fillRandom, contains, setDifference, sort, flip, readFirstFrame, readNextFrame, loadPosition, loadPositionCenter, ...** 

## How to compile the code

```bash
cd Analysis
./configure path/to/<sourceCode>.cpp  # compiler one of the C++ source file <courceCode>.cpp
```
