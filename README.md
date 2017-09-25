# BnB_CETSP

A branch-and-bound algorithm for the CETSP.

This is a prototype code that has not been developed for distribution.

## Building an executable

To build an executable, run the following command from the project root directory "BnB_CETSP":
```
make
```
Please, note that:
 * You will need to edit the makefile in order to fix the correct path to your CPLEX installation.
 * You will need to install the "GMP MP Bignum Library" (https://gmplib.org/) in order to run this code.

## Calling command:

The code can be run by using the following command structure
```
./exeCETSP [Path to the instance] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE]
```
Where:

 * OPTIONS = "2D" or "3D"
 * OVERLAP FACTOR = One of the following values "{0.1, 0.25, 0.5, 1.0, 1.5}"
 * TIME LIMIT = Any integer value
 * BRANCHING RULE = One of the following values "{V1,SB}"
 * BRANCHING Strategy = One of the following values "{DFS,BFS,BeFS}"
 * ROOT SELECTION = ONe of the following values "{1,2,3}"
 * S.B SIZE = Any integer value

Not all combinations of parameters will work for all instances. Check the paper for details (https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2016.0711).

Examples:

Running a 2D instance:
```
./exeCVXHULL 2D/rotatingDiamonds1.txt 2D 1.0 3600 V1 BeFS 1 1
```

Running a 3D instance:
```
./exeCVXHULL 3D/d493.txt 3D 0.5 3600 V1 BeFS 1 1
```

## Authors

* **Walton Pereira Coutinho** - *Initial work* - [BnBCETSP](https://github.com/waltonpcoutinho/BnB_CETSP)

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details


