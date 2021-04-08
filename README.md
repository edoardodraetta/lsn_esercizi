# lsn_esercizi
LSN exercise delivery.

#### RAN NYU

La routine é RANdom New York University (RANNYU). Scritta negli anni 80 con un moltiplicatore suggerito nel 2 volume di Knuth "Art of scientific computing" p 102. Riadattata poi a inizio anni 90 nel gruppo di M.H. Kalos (da parte del Dr. Panoff). 

Si tratta di un Linear Congruential Generator with prime addend a 48 bit. 
La relazione che genera x_{n+1} da x_n é la seguente: x_{n+1}= ( a x_{n} + c ) mod m 
dove a é il moltiplicatore, c é l'addendo (numero primo con certe proprietà), m é il modulo, in questo caso m=2\*\*48. 
Con questo moltiplicatore e questo modulo il periodo del generatore risulta P=2\*\*48~10\*\*14. Gli x_n sono numeri interi e corrispondono allo stato salvato nella routine. Il numero random reale in [0,1) si ottiene con x_n/m. 

Notare che per ogni c si ha un generatore diverso. In un lavoro di Percus e Kalos (che vi ho aggiunto) sono state mostrate le proprietà che due diversi c (c1 e c2) devono avere affinché le serie generate (x1_n e x2_n) siano effettivamente indipendenti. I 384 numeri primi che già avete in Primes soddisfano questo criterio. Seguendo quel criterio, un mio collaboratore (Dr. G. Bertaina) me ne ha ottenuti 32001. Li trovate in primes32001.in, in ordine decrescente. 

Notate che nella routine lo stato a 48 bit é salvato in 4 interi minori di 4096 (=2\*\*12). Il tutto, comprese le operazioni, é scritto in base 4096 (anche in numeri primi nel file Primes, per questo sono coppie)

### Exercise 1

Run `make` in `./rng/` to generate `main.exe`, which will run all experiments in exercise 1. 

The jupyter notebook needs to be run sequentially. 

#### TO DO:

 - [ ] Is the Cauchy fit a true fit?
 - [x] chi2
 - [ ] Double check python code
 - [ ] Pass average calculation by reference?
 - [x] move data folder out of rng
 - [ ] 2d theta correct

### Exercise 2

 - [x] Importance Sampling
 - [ ] Are my RWs correct?
 - [x] RW fit

### Exercise 3

 - [x] do the math/scratchwork before you code!
 - [ ] implement max comparator
 - [ ] reinitialize rng
 - [ ] implement call/puts (same file!)  
