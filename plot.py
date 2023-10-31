import numpy as np
import matplotlib.pyplot as plt

with open('toPlot.txt', 'r') as file:
    content = file.read()
    content = content.split('\n')
    content = [x.split(' ') for x in content]
    x = [float(i[0]) for i in content if i[0] != '']
    y = [float(i[1]) for i in content if i[0] != '']
    z = [float(j[-1]) for j in content if j[-1] != '']

plt.scatter(x, y, marker='o', s=0.03)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
#plt.zliabel("z [m]")
plt.title("Ruch cząstki")
plt.savefig("ruch.png")
plt.close()

with open('Potencjaly', 'r') as file:
    content = file.read()
    content = content.split('\n')
    content = [x.split(' ') for x in content]

content = [[float(j) for j in i if j != ''] for i in content if i != ['']]
x = np.arange(1, len(content[0])+1)
y = np.arange(1, len(content)+1)
X, Y = np.meshgrid(x,y)
#print(f"len(X) = {len(X)}, lenY = {len(Y)}, len(content[0]) = {len(content[0])}, len(content) = {len(content)}")

plt.contourf(X, Y, content, cmap='cool')
plt.title('Potencjał dla indeksu Z=1, po 1000 interacji.')
plt.xlabel('oś x [kratka]')
plt.ylabel('oś y [kratka]')
plt.colorbar()
plt.savefig("pot.png")
plt.close()
