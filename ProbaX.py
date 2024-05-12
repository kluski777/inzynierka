#!/usr/bin/env python
# coding: utf-8

# In[27]:


# obróbka danych
import numpy as np
import matplotlib.pyplot as plt

path = "C:\\Users\\klusk\\cpp\\quasiGlupoty\\"
# path2 = "C:\\Users\\klusk\\cpp\\quasiGlupoty\\"
path2 = "D:\\MOFIT\\inzynierka\\New folder\\"

def read4D(fileName):
    with open(fileName, "rb") as content:
        content = content.read().decode('utf-8')
        content = content.strip().split('\r')
        content = [i.strip().split('\t') for i in content]
        content = [[j.strip().split('\n') for j in i] for i in content]
        content = [[[k.strip().split(' ') for k in j] for j in i] for i in content]

    content = [[[[float(l) for l in k] for k in j] for j in i] for i in content]
    print(f"Wymiary 4D {fileName}: {len(content)} x {len(content[0])} x {len(content[0][0])} x {len(content[0][0][0])}")
    return np.array(content)

def read3D(fileName):
    with open(fileName, "rb") as content:
        content = content.read().decode('utf-8')
        print("Jestem 1")
        content = content.strip().split('\t')
        print("Jestem 2")
        content = content[:2000]
        print("Jestem 3")
        content = [i.strip().split('\n') for i in content]
        print("Jestem 4")
        content = [[j.strip().split(' ') for j in i] for i in content]
    print("Jeszcze chwilka")
    content = [[[float(k) for k in j if k != ''] for j in i] for i in content]
    print(f"Wymiary 3D {fileName}: {len(content)} x {len(content[0])} x {len(content[0][0])}")
    return np.array(content)

def read2D(fileName):
    with open(fileName, "rb") as file:
        content = file.read().decode('utf-8')
        content = content.strip().split('\n')
        content = [i.strip().split(' ') for i in content]

    content = [[float(j) for j in i] for i in content]
    print(f"Wymiary 2D {fileName}: {len(content)} x {len(content[0])}")
    return np.array(content)


# In[303]:


# ### NEWTON
# dt, Nx, Nz, dx, itersy, coIleRest, coIlePot = read2D('info.bin')[0]
# czas = np.arange(itersy/coIleRest)*(coIleRest*dt/60/60/24/365)
# Nx, Nz, itersy, coIleRest, coIlePot= int(Nx), int(Nz), int(itersy), int(coIleRest), int(coIlePot)

# potN = read4D(f'{path}Potencjal.bin')
# # den = read4D("gestosc")
# ra  = read3D(f"{path}Czatki.bin") # to są te z Newtona liczone
# En  = read2D(f"{path}Energia.bin")

# rxN  = ra[:, :, 0]
# ryN  = ra[:, :, 1]
# axN  = ra[:, :, 2]
# ayN  = ra[:, :, 3]
# ekN  = En[:, 0]
# epN  = En[:, 1]
# LN   = En[:, 2]


# In[28]:


### GAUSS
dt, Nx, Nz, dx, itersy, coIleRest, coIlePot, N = read2D(f'{path2}info')[0]
czas = np.arange(itersy/coIleRest)*(coIleRest*dt/60/60/24/365)
T = czas[:]/2.2e8
Nx, Nz, itersy, coIleRest, coIlePot, N = int(Nx), int(Nz), int(itersy), int(coIleRest), int(coIlePot), int(N)
nx, nz = (Nx-1)/2, (Nz-1)/2

pot = read4D(f"{path2}Potencjal")
den = read4D(f"{path2}gestosc")
ra  = read3D(f"{path2}Czatki")
En  = read2D(f"{path2}energie")

rx  = ra[:, :, 0]
ry  = ra[:, :, 1]
ax  = ra[:, :, 2]
ay  = ra[:, :, 3]
ek  = En[:, 0]
ep  = En[:, 1]
L   = En[:, 2]

czas = np.arange(len(rx[:,:]))*(coIleRest*dt/60/60/24/365)
ile = len(ek)
ro = 6e42/int(np.sum(den[0, :, :, :]))/dx**3


# In[259]:


len(rx[:,:])


# In[29]:


save2 = "D:\\MOFIT\\inzynierka\\Więcej Cząstek\\Okrąg\\"

# X, Z = np.meshgrid(z, x)

z = np.arange(Nz) - (Nz-1)/2

X, Y = np.meshgrid((np.arange(Nx) - (Nx-1)/2)*dx, (np.arange(Nx) - (Nx-1)/2)*dx)

for i in range(int(Nz)):
    plt.title(f"Dla płaszczyzny z = {i} na końcu symulacji")
    plt.contourf(Y, X, pot[0, :, :, i], levels=100)
    plt.ylabel('Położenie [m]')
    plt.xlabel('Położenie [m]')
    plt.colorbar()
    plt.show()


# In[30]:


from matplotlib.colors import LinearSegmentedColormap
colors = [(1, 1, 1), (0.0, 0.0, 0.0)]  # Light gray to black
cmap = LinearSegmentedColormap.from_list('custom_gray', colors, N=256)

for j in range(int(itersy/coIlePot)):
    for i in range(1, int(Nz)-1):
        plt.title(f'Gęstość masy po {czas[int(j*coIlePot/coIleRest)]:.2e} latach')
        plt.contourf(Y, X, den[j, :, :, i]*ro, cmap=cmap, levels=50)
        plt.colorbar()
        plt.xlabel('Położenie [m]')
        plt.ylabel('Położenie [m]')
        plt.savefig(save2+'ro_z'+str(i)+'_iter'+str(j))
        plt.pause(.1)


# In[31]:


# narysować siatke jeszcze.
pion   = np.zeros([Nx-1, Nx-1], dtype=np.float64)
poziom = np.zeros([Nx-1, Nx-1], dtype=np.float64)
x      = (np.arange(Nx-1) - Nx/2)*dx + dx
for i in range(Nx-1):
    for j in range(Nx-1):
        pion[i, j] = x[i]
        poziom[i, j] = x[i]


# In[32]:


plt.title('Energia kinetyczna, potencjalna i całkowita')
plt.plot(czas, ek, label=r'$T$')
plt.plot(czas, ep, label=r'$V$')
plt.plot(czas, ek + ep, label=r'$E = T + V$')

plt.ylabel('Energia [J]')
plt.xlabel('Czas [lata]')
plt.legend(loc='lower right')
plt.savefig(save2+'energie')
# cząstki, co wypadły nie wpływają na potencjał


# In[33]:


plt.plot(czas, L, label='metodą PiC') # weź to dokładniej walnij, może na Eulerze
# plt.plot(czas, LN, label='metodą DNB')
plt.title("Moment pędu")

plt.ylim(2.4e67, 3.5e67)
plt.ylabel(r'Moment pędu [$\frac{kg m^2}{s^2}$]')
plt.xlabel('Czas [lata]')
plt.savefig(save2+'L')
# plt.legend()


# In[ ]:


int(np.sum(den))


# In[ ]:





# In[34]:


plt.plot(czas, ek, label='metodą PiC')
# plt.plot(czas, ekN, label='metodą DNB')
plt.title("Energia kinetyczna")
plt.ylabel('E [J]')
plt.xlabel('Czas [lata]')
plt.savefig(save2+'ek')


# In[35]:


plt.title("Energia całkowita")
plt.plot(czas, ek + ep, label='metodą PiC')

plt.xlabel('Czas [lata]')
plt.ylabel('Energia [J]')
plt.savefig(save2+'Ec')


# In[ ]:


len(rx[:,0])


# In[36]:


cnt = 0
since = 110
to = 200
plt.title(f'Położenia cząstek w czasie od {czas[since]:.2e} do {czas[to-1]:.2e}')
for i in range(N):
    if np.sum(np.where(np.abs(rx[since:to, i]) < nx*dx, 0, 1)) == 0 and np.sum(np.where(np.abs(ry[since:to, i]) < nx*dx, 0, 1)) == 0:
        plt.plot(rx[since:to, i], ry[since:to, i], label='Czastka '+str(i+1)+' PiC', linewidth=.1)
    else:
        cnt += 1
        
# for i in range(N):
#     plt.plot(rxN[:, i], ryN[:, i], label='Czastka '+str(i+1)+' DNB', linewidth=.5)

for i in range(Nx-1):
    plt.plot(pion[i, :], x, linewidth=.1, color='b')
    plt.plot(x, poziom[i, :], linewidth=.1, color='b')
    
plt.xlim([np.min(x), np.amax(x)])
plt.ylim([np.min(x), np.amax(x)])
# plt.legend(loc='upper right')
plt.ylabel('x [m]')
plt.xlabel('y [m]')
plt.axis('equal')
plt.savefig(save2+'ruchy')

print(f"{int(cnt/N*100)}%")


# In[ ]:


ax[:]


# In[37]:


for i in range(N):
    plt.title('Położenie we współrzędnej x.')
#     plt.plot(czas, rxN[:, i], label='Czastka Gauss ' + str(i), linewidth=1)
    plt.plot(czas, rx[:, i], label='Czastka Gauss ' + str(i), linewidth=1)
    
plt.xlabel('czas [lata]')
plt.ylabel('Wychyelenie [m]')
plt.savefig(save2+'x')
# plt.legend(loc='upper right')


# In[38]:


plt.title('Położenie we współrzędnej y.')
for i in range(N):
    if np.sum(np.where(np.abs(ry[:, i]) < 5e21, 0, 1)) == 0:
        plt.plot(czas, -ry[:, i], label='Czastka '+str(i+1)+' PiC', linewidth=1)

# for i in range(N):
#     plt.plot(czas, ryN[:, i], label='Czastka '+str(i+1)+' DNB', linewidth=1)

# plt.legend(loc="upper right")
plt.xlabel('czas [lata]')
plt.ylabel(r'Wychyelenie w kierunku $y$ [m]')
plt.savefig(save2+'y')


# In[39]:


plt.title(r'$a_x$, dla ilości kratek w każdym z kierunków' + f'\n(Nx, Ny, Nz) = ({Nx}, {Nx}, {Nz})' + f' dt = {dt:.1e}')
for i in range(N):
    plt.plot(czas, ax[:, i], label='Czastka '+str(i+1)+' PiC', linewidth=.75)
# for i in range(N):
#     plt.plot(czas, axN[:, i], label='Czastka '+str(i+1)+' DNB', linewidth=.75)
# plt.legend(loc="lower right")
plt.ylabel(r'$a_x$ [$\frac{m}{s^2}$]')
plt.xlabel('Czas [lata]')
plt.savefig(save2+'ax')


# In[40]:


plt.title(r'$a_y$, dla ilości kratek w każdym z kierunków' + f'\n(Nx, Ny, Nz) = ({Nx}, {Nx}, {Nz})' + f' dt = {dt:.1e}')
for i in range(N):
    plt.plot(czas, ay[:, i], label='Czastka '+str(i)+' PiC', linewidth=0.75)
# for i in range(N):
#     plt.plot(czas, ayN[:, i], label='Czastka '+str(i)+' DNB', linewidth=0.75)
    
plt.ylabel(r'$a_y$ [$\frac{m}{s^2}$]')
plt.xlabel('Czas [lata]')
plt.savefig(save2+'ay')
# plt.legend(loc='lower left')


# In[41]:


plt.title('Kwadrat przyspieszenia.')
for i in range(N):
    plt.plot(czas, np.sqrt(ay[:, i]**2+ax[:, i]**2), label='Czastka Gauss '+str(i), linewidth=1)
# plt.legend(loc="lower left")
plt.ylabel(r'$\sqrt{\vec{a}^2}$ [$\frac{m}{s^2}$]')
plt.xlabel('Czas [lata]')
plt.savefig(save2+'a2')


# In[ ]:


np.sqrt(ayN[:, i]**2+axN[:, i]**2)


# In[42]:


plt.title(r'$v_x$')
for i in range(N):
    plt.plot(czas[:-1], np.diff(rx[:, i])/dt, label='Czastka '+str(i), linewidth=.75) # chyba zbyt często są zbierane wyniki
plt.xlabel('Czas [lata]')
plt.ylabel(r'$v_x$ [$\frac{m}{s}$]')
plt.savefig(save2+'vx')
# plt.legend(loc='lower left')


# In[ ]:





# In[43]:


plt.title(r'$v_y$')
for i in range(N):
    plt.plot(czas[:-1], np.diff(ry[:, i])/dt, label='Czastka '+str(i), linewidth=.75)
plt.ylabel(r'$v_y$ [$\frac{m}{s}$]')
plt.xlabel('czas [lata]')
plt.savefig(save2+'vy')
# plt.legend(loc='lower left')


# In[44]:


plt.title(r'$\vec{v}^2$')
for i in range(N):
    plt.plot(czas[:-1], np.sqrt(np.diff(ry[:, i])**2 + np.diff(rx[:, i])**2)/dt, label='Czastka '+str(i+1)+' PiC', linewidth=0.5)

plt.xlabel('time [yrs]')
# plt.legend(loc='lower left')
plt.ylabel(r'v [$\frac{m}{s}$]')
plt.savefig(save2+'v2')


# In[45]:


meanR = np.zeros([int(itersy/coIleRest)])
median = np.zeros([int(itersy/coIleRest)])
stddevR = np.zeros([int(itersy/coIleRest)])
turet = np.zeros(N)
for i in range(int(itersy/coIleRest)):
    for j in range(N):
        turet[j] = np.sqrt(rx[i, j]**2 + ry[i, j]**2)
    meanR[i] = np.mean(turet)
    median[i] = np.median(turet)
    stddevR[i] = np.std(turet)

plt.title('Średnie r')
plt.plot(czas, meanR)

plt.xlabel('time [yrs]')
plt.ylabel('Odległości od środka [m]')
plt.savefig(save2+'avgr')


# In[ ]:





# In[46]:


ileJestW = np.zeros(int(itersy/coIleRest), dtype=np.float64)
for i in range(int(itersy/coIleRest)):
    ileJestW[i] = np.sum(np.where((np.abs(rx[i, :]) < nx*dx) & (np.abs(ry[i,:]) < nx*dx), 0, 1))*100/N

fig, ax1 = plt.subplots()
ax1.set_title('mediana odległości i część cząstek które wyszły za siatkę')
ax1.plot(czas, median, label='mediana')
ax1.set_xlabel('time [lata]')
ax1.set_ylabel('Odległości od środka [m]')

ax2 = ax1.twinx()
ax2.plot(czas, ileJestW, color='red', label='Cząstki które wypadły')
ax2.set_ylabel('Część cząstek poza siatką [%]')
ax2.legend(loc='upper left')
ax1.legend(loc='lower right')
plt.savefig(save2+'medianR')


# In[47]:


meanV = np.zeros([int(itersy/coIleRest)-1])
stddevV = np.zeros([int(itersy/coIleRest)-1])
turet = np.zeros(N-1)
for i in range(int(itersy/coIleRest)-1):
    for j in range(N-1):
        turet[j] = np.sqrt((rx[i+1, j]-rx[i, j])**2 + (ry[i, j] - ry[i+1, j])**2)/dt
    stddevV[i] = np.std(turet)
    meanV[i] = np.mean(turet)

plt.title(r'Średnie v')
plt.plot(czas[:-1], meanV)

plt.xlabel('time [yrs]')
plt.ylabel(r'$\sqrt{\vec{v}^2} [\frac{m}{s}]$')

plt.savefig(save2+'avgv')


# In[48]:


plt.title(r'$\sigma_v$')
plt.plot(czas[:-1], stddevV)

plt.xlabel('time [yrs]')
plt.ylabel(r'$\sigma(v)     [\frac{m}{s}]$')
plt.savefig(save2+'sigmav')


# In[49]:


meana = np.zeros([int(itersy/coIleRest)])
stda = np.zeros([int(itersy/coIleRest)])
turet = np.zeros(N)
for i in range(int(itersy/coIleRest)):
    for j in range(N):
        turet[j] = np.sqrt(ax[i, j]**2 + ay[i, j]**2)
    meana[i] = np.median(turet)
    stda[i] = np.std(turet)

plt.title(r'Średnie $\vec{a}$, dla (Nx, Nx, Nz) = ' + f'({Nx},{Nx},{Nz}) i dt = {dt:.0} s')
plt.plot(czas, meana)

plt.xlabel('time [yrs]')
plt.ylabel(r'$\sqrt{\vec{a}^2} [\frac{m}{s}]$')
plt.savefig(save2+'avga')


# plt.title(r'$\sigma(\vec{a})$, dla (Nx, Nx, Nz) = ' + f'({Nx},{Nx},{Nz}) i dt = {dt:.0} s')
# plt.plot(czas, stda)
# 
# plt.xlabel('time [yrs]')
# plt.ylabel(r'$\sqrt{\vec{a}^2} [\frac{m}{s}]$')
# plt.savefig(save2+'stda')
