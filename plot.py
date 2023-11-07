import numpy as np
import matplotlib.pyplot as plt

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
    print(f"len(X) = {len(X)}, lenY = {len(Y)}, len(content[0]) = {len(content[0])}, len(content) = {len(content)}")

    plt.contourf(X, Y, content, cmap='cool', levels=100)
    plt.title(pltTitle)
    plt.xlabel('oś x [kratka]')
    plt.ylabel('oś y [kratka]')
    plt.colorbar()
    plt.savefig(fileName)
    plt.close()

def plotPath(fileName, pltTitle):
    with open(fileName, "r") as content:
        content = content.read()
        content = content.split('\n')
        content = [i.split(' ') for i in content if i != '']
 
    content = [[float(j) for j in i if j != ''] for i in content]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    content = np.array(content)

    x = [i[0] for i in content]
    y = [i[1] for i in content]
    z = [i[2] for i in content]

    ax.plot(x, y, z)
    ax.set_title(pltTitle)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    fig.savefig(fileName + '3d')
    plt.close(fig)


    plt.plot(x, y)
    plt.title(pltTitle)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.savefig(fileName + '2d')
    plt.close()

#plotPot("Potencjal1", "Potencjał na całej siatce dla Z_indx = 0")
#plotPot("Potencjal2", "Potencjał na całej siatce dla Z_indx = 1")
#plotPot("Potencjal3", "Potencjał na całej siatce dla Z_indx = 2")
#plotPot("Potencjal4", "Potencjał na całej siatce dla Z_indx = 3")

plotPath("czastka1", "Droga cząstki 1")
plotPath("czastka2", "Droga cząstki 2")
plotPath("czastka3", "Droga cząstki 3")
plotPath("czastka4", "Droga cząstki 4")
