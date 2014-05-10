# DEFINIRE VARIABILI, MATRICI, PARAMETRI
u = [[0]*dim]*dim

# CONDIZIONI A CONTORNO DEL SISTEMA


# FORMULA Gauss-Jacobi

u_new[i][j] = (u_old[i-1][j] + u_old[i+1][j] + u_old[i][j-1] + u_old[i][j-1] \
               + h**2*roh[i][j])*1./4.

u_old = u_new


# FORMULA Gauss-Seidel

u[i][j] = 1./4.*(u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j-1] + h**2*roh[i][j])


# FORUMULA S.O.R.

u[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j-1] + h**2*roh[i][j])\
          *omega/4. + (1 - omega)*u[i][j]


# DETERMINAZIONE ITERAZIONI NECESSARIE ALLA CONVERGENZA


# DETERMINAZIONE OMEGA OTTIMALE METODO S.O.R.
