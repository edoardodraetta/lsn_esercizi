
### Exercise 9
1. Perform a random search with 10, 15, 32 cities on a circle
2. Perform a genetic search with 10, 15, 32 cities on a circle
3. Perform a random search with 10, 15, 32 cities on a square
4. Perform a g search with 10, 15, 32 cities on a square

#### Data structures

- Each path is a row vector of integers `urowvec`
- The population will be represented by a `field<urowvec>`, such that armadillo methods will be able to called on each path
- Each path will be assumed to start and end at the first city
