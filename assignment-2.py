# DEFINIRE VARIABILI, MATRICI, PARAMETRI
dim = 100
u = [[0 for i in range(dim)] for i in range(dim)]
rho = [[0 for i in range(dim)] for i in range(dim)]
#rho[int(0.8*dim)][int(0.8*dim)] = 5
#rho[int(0.2*dim)][int(0.2*dim)] = -5
file_jacobi = open('jacobi.dat', 'w')
file_gauss = open('gauss.day', 'w')
file_SOR = open('SOR.dat', 'w')
# CONDIZIONI A CONTORNO DEL SISTEMA
for pos in range(dim):
    u[0][pos] = 100
    u[dim-1][pos] = 100
    u[dim-1][0] = -150
    u[pos][dim-1] = -500


u_new = u_old = u
h = 0.1
# FORMULA Gauss-Jacobi
converg = False
sum_new = sum_old = 0
diff_old = 0
step = 0
while not converg:
    for i in range(1,dim-1):
        for j in range(1,dim-1):
            u_new[i][j] = (u_old[i-1][j] + u_old[i+1][j] + u_old[i][j-1] \
                       + u_old[i][j+1] + h**2*rho[i][j])*1./4.
            sum_new += u[i][j]
            
    diff_new = sum_new - sum_old
    converg = abs(diff_new - diff_old) == 0
    step += 1
    print diff_new-diff_old, "\t", sum_old, "\t", sum_new, "\t", converg, "\t#", step
    u_old = u_new
    sum_old = sum_new
    diff_old = diff_new

for i in range(dim):
    for j in range(dim):
        file_jacobi.write(str(i))
        file_jacobi.write("\t")
        file_jacobi.write(str(j))
        file_jacobi.write("\t")
        file_jacobi.write(str(u_new[i][j]))
        file_jacobi.write("\n")


a = raw_input(">>>")
# DEFINIRE VARIABILI, MATRICI, PARAMETRI
dim = 50
u = [[0 for i in range(dim)] for i in range(dim)]
rho = [[0 for i in range(dim)] for i in range(dim)]
#rho[int(0.8*dim)][int(0.8*dim)] = 5
#rho[int(0.2*dim)][int(0.2*dim)] = -5
file_jacobi = open('jacobi.dat', 'w')
file_gauss = open('gauss.day', 'w')
file_SOR = open('SOR.dat', 'w')
# CONDIZIONI A CONTORNO DEL SISTEMA
for pos in range(dim):
    u[0][pos] = 100
    u[dim-1][pos] = 100
    u[dim-1][0] = -150
    u[pos][dim-1] = -500


h = 0.1
# FORMULA Gauss-Jacobi
converg = False
sum_new = sum_old = 0
diff_old = 0
step = 0
# FORMULA Gauss-Seidel

while not converg:
    for i in range(1,dim-1):
        for j in range(1,dim-1):
            u[i][j] = 1./4.*(u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] \
                             + h**2*rho[i][j])
            sum_new += u[i][j]
            
    diff_new = sum_new - sum_old
    converg = abs(diff_new - diff_old) == 0
    step += 1
    print diff_new-diff_old, "\t", sum_old, "\t", sum_new, "\t", converg, "\t#", step
    sum_old = sum_new
    diff_old = diff_new

for i in range(dim):
    for j in range(dim):
        file_gauss.write(str(i))
        file_gauss.write("\t")
        file_gauss.write(str(j))
        file_gauss.write("\t")
        file_gauss.write(str(u[i][j]))
        file_gauss.write("\n")


a = raw_input(">>>")
# DEFINIRE VARIABILI, MATRICI, PARAMETRI
dim = 50
u = [[0 for i in range(dim)] for i in range(dim)]
rho = [[0 for i in range(dim)] for i in range(dim)]
#rho[int(0.8*dim)][int(0.8*dim)] = 5
#rho[int(0.2*dim)][int(0.2*dim)] = -5
file_jacobi = open('jacobi.dat', 'w')
file_gauss = open('gauss.day', 'w')
file_SOR = open('SOR.dat', 'w')
# CONDIZIONI A CONTORNO DEL SISTEMA
for pos in range(dim):
    u[0][pos] = 100
    u[dim-1][pos] = 100
    u[dim-1][0] = -150
    u[pos][dim-1] = -500


h = 0.1
# FORMULA Gauss-Jacobi
converg = False
sum_new = sum_old = 0
diff_old = 0
step = 0
omega = 1.91
# FORUMULA S.O.R.

while not converg:
    for i in range(1,dim-1):
        for j in range(1,dim-1):
            u[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] \
                       + h**2*rho[i][j])*omega/4. + (1 - omega)*u[i][j]
            sum_new += u[i][j]
            
    diff_new = sum_new - sum_old
    converg = abs(diff_new - diff_old) == 0
    step += 1
    print diff_new-diff_old, "\t", sum_old, "\t", sum_new, "\t", converg, "\t#", step
    sum_old = sum_new
    diff_old = diff_new

for i in range(dim):
    for j in range(dim):
        file_SOR.write(str(i))
        file_SOR.write("\t")
        file_SOR.write(str(j))
        file_SOR.write("\t")
        file_SOR.write(str(u[i][j]))
        file_SOR.write("\n")


# DETERMINAZIONE ITERAZIONI NECESSARIE ALLA CONVERGENZA


# DETERMINAZIONE OMEGA OTTIMALE METODO S.O.R.


file_jacobi.close()
file_gauss.close()
file_SOR.close()

