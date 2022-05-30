import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import NonLinear.chaos as chaos


uGrid,wGrid = np.linspace(-1,1,400),np.linspace(-1,1,400)
eta = 0.92
itrMap = -1*np.ones((400,400))

for u_i,u in enumerate(uGrid):
    for w_i,w in enumerate(wGrid): 
        if ( u ** 2 + w **2 > 1):
            # print(u,w,u**2 + w**2)
            itrMap[w_i][u_i] = 0
            continue
        
        itr = 0
        u_prev,w_prev = 0,0
        # u,w = 0.5,0.5
        U,W = u,w
        while(itr < 99):
            U_next,v,W_next = chaos.forward(U,0,W,eta)
            itr += 1
            
            if itr % 2 :
                if ( abs(U_next - u_prev) < 0.005 and abs(W_next- w_prev) < 0.005):
                    if ( abs(U) < 0.05 and abs(W) < 0.05 ):
                        itrMap[w_i][u_i] = 1
                    elif ( U - W > 0 and U > 0.5):
                        itrMap[w_i][u_i] = 2
                    elif (W - U > 0 and W > 0.5):
                        itrMap[w_i][u_i] = 3
                    else:
                        itrMap[w_i][u_i] = 1
                    break
            u_prev,w_prev = U,W
            U,W = U_next,W_next


          

colors = ['white', 'red','blue','lightblue']
bounds = [-0.5,0.5,1.5,2.5,3.5]

cmap = mpl.colors.ListedColormap(colors)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

plt.imshow(itrMap, origin='lower', interpolation='none', cmap=cmap, norm=norm)
plt.show()
plt.savefig('asd.png')
