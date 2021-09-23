
The two steps to obtain syntenic information for each intact LTR-RT:

First, obtain pairwise syntenic information for each LTR-RT between any two genomes.
	- There are 27*26/2 = 351 pairs.
	- Each pair took ~ one month of a server (36 cpu, 360G memmory) and generated approx. 2T of intermediate data.
	- The final results are small.

Second, combine pairwise information into a pan-genome syntenic matrix.
	- convert approximate matches into exact intact LTR-RT coordinates
	- combine each pairwise genome pairs
	- resolve overlaps and conflicts
	- This step took ~ one day.
