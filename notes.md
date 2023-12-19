# General
- You lack a readme in which you explain how to compile the code and what to expect
- The structure of the repo is a bit poor, it would have been better to have folders with source files, tests, scripts and so on.
- Use a `.gitignore` to ignore generated files 
- Nice that you have many tests, the interface of the Multigrid and Lattice classes is well done

# Code
## Major
- The Makefile rules should separate the compilation from the execution
- it is not a good practice to pollute the global scope with `using namespace std;`
- if a class is not template try to split into header and souce files
- check the other notes in the code

## Minor
- The octave script could save the report on a file