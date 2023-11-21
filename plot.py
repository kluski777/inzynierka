import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def plotPot(fileName, pltTitle):
    with open(fileName, 'r') as file:
        content = file.read()
        content = content.split('\n')
        content = [x.split(' ') for x in content if x != [''] or x != '']
    content = content[:-1]

    content = [[float(j) for j in i if j != ''] for i in content if i != [''] or i != '' or i != []]
    x = np.arange(1, len(content[0])+1)
    y = np.arange(1, len(content)+1)
    X, Y = np.meshgrid(x,y)
    print("Potencjal")
    #print(f"len(X) = {len(X)}, lenY = {len(Y)}, len(content[0]) = {len(content[0])}, len(content) = {len(content)}")

    plt.contourf(X, Y, content, cmap='cool', levels=100)
    plt.title(pltTitle)
    plt.xlabel('oś x [kratka]')
    plt.ylabel('oś y [kratka]')
    plt.colorbar()
    plt.savefig(fileName)
    plt.close()

def read3D(fileName):
    with open(fileName, "rb") as content:
        content = content.read().decode('utf-8')
        print(f"len(content) = {len(content)}")
        content = content.split('\t')
        content = content[:-1]
        content = [i.strip().split('\n') for i in content]
        content = [[j.strip().split(' ') for j in i] for i in content]
        
    content = [[[float(k) for k in j] for j in i] for i in content]

    return content

plotPot("Potencjal1", "Potencjał na całej siatce dla Z_indx = 0")
plotPot("Potencjal2", "Potencjał na całej siatce dla Z_indx = 1")
plotPot("Potencjal3", "Potencjał na całej siatce dla Z_indx = 2")
plotPot("Potencjal4", "Potencjał na całej siatce dla Z_indx = 3")
plotPot("densityBef", "Gęstość dla Z_indx=1 przed iteracją")
plotPot("densityAfter", "Gęstość dla Z_indx=1 po iteracji")


wholeData = read3D("droga")
wholeData = np.array(wholeData)
toPlot = wholeData[0]

fig, ax = plt.subplots()

scatter = ax.scatter(toPlot[:,0], toPlot[:,1], s=1)


def update(frame):
    toPlot = wholeData[frame]
    scatter.set_offsets(np.column_stack((toPlot[:,0], toPlot[:,1])))

    return scatter,

animation = FuncAnimation(plt.gcf(), update, frames=40, interval=1, blit=True, repeat=False)

animation.save('scatterAnimation.gif', writer='pillow')

