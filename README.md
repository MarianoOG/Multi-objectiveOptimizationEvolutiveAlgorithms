# Multi-objective Optimization Evolutive Algorithms

A wide variety of real world problems have several (often conflicting) objectives that need to be optimized at the same time. They are called multi-objective optimization problems (MOPs) and their solution involves finding a set of decision variables that represent the best trade-offs among all the objectives.

In this project, I programmed NGSA-III algorithm to solve multi-objective problems with constrains. I also programmed some predefined problems from the DTLZ set in order to compare the results against some other authors and algorithms. The results were close to the original paper, even though some of the parts of the algorithm were hidden by the author.

## Applications

The real world problem application was the design of an electric circuit with objectives were the following: maximize: f(1) - Open Loop Gain, f(2) - Unitary Gain Frequency, f(3) - Common Mode Rejection Rate, f(4) - Slew Rate; minimize: f(5) - Input Current, f(6) - Area. The results were sub optimal, because NSGA-II generated better outcomes for this problem, this is probably because of the lack of the mentioned sections.

## More information

These are the significate files in which I contributed during a summer internship done in 2016 at the Computer Science Department of CINVESTAV at the *Verano de la Investigación Científica y Tecnológica del Pacífico XXI*. Assessors: PhD. Carlos A. Coello Coello (ccoello@cs.cinvestav.mx) and PhD. Candidate Raquel Hernández Gómez.

Please refer to the [presentation](MOEA.pdf) for more information.
