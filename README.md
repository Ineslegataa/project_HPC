High Performance Computing – Exam 2024

Questa repository contiene gli assignment svolti per il corso di High Performance Computing, parte del MSc in Data Science and Scientific Computing presso SISSA (2024).
Il lavoro è incentrato sull’implementazione di algoritmi scientifici paralleli e sull’analisi delle prestazioni tramite MPI, OpenMP e approcci ibridi.
Contenuto della repository

assignment_1/

    Operazione stencil a cinque punti in 2D

    Implementata con OpenMP (memoria condivisa)

    Benchmark su diversi numeri di thread

    Analisi di speedup ed efficienza

assignment_2/

    Esercizio 2c: generazione dell’insieme di Mandelbrot

        Implementato in C++

        Parametrico da linea di comando:

        ./mandelbrot n_x n_y x_L y_L x_R y_R I_max

        Output generato in formato .pgm in scala di grigi (mandelbrot.pgm)

        Pronto per estensioni OpenMP/MPI

    Grafici di latenza (naive_model.png) per diversi algoritmi MPI collettivi:

        basic_linear, pipeline, binomial, scatter_allgather, ecc.

        Analisi della latenza (in microsecondi) al variare del numero di core

Compilazione e utilizzo

Per compilare:

make

Per eseguire:

./mandelbrot n_x n_y x_L y_L x_R y_R I_max

Note

Tutti i contenuti sono stati sviluppati individualmente per l’esame di HPC 2023/2024.
Sono condivisi solo a scopo accademico. Il riutilizzo non autorizzato è scoraggiato.
