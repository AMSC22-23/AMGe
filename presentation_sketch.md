Introduzione
Quale problema si vuole risolvere (laplaciano, condizioni di dirichlet al bordo)
Quali semplificazioni sono state fatte (mesh quadrata, numero di suddivisioni)
Essenza del metodo multigrid
Cosa non è multigrid (le due false partenze)

Come abbiamo inizialmente usato le feature del c++ per generalizzare il concetto di mesh (Mesh.hpp)
Riferimento all'inner nodes, e all'indicizzazione

Di conseguenza come i metodi di calcolo sono agnostici rispetto alla mesh, devono solo sapere quali sono i nodi interni e dove sono i vicini

Presentazione dei due smoother jacobi e gauss seidel, descrivere particolarità soprattutto riferendosi alla facilità di parallelizzazione

Mostrare il benchmark che giustifica il fatto di usare openmp (limite teorico di speedup raggiunto)

Parlare del bug che ci ha tenuto fermi due settimane, di come si sia prodotto dal fatto di lavorare matrix-free e dal fatto di voler fare le cose sofisticate (divisioni diverse nella x e nella y)

Descrivere il tradeoff nel multigrid a più livelli tra convergenza nel minor numero di iterazioni e velocità di calcolo delle singole iterazioni, corredato da benchmark
