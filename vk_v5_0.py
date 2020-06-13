#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gamma_function
import math


# In[ ]:





# In[2]:


class Cuckoo():
    
    X_min = 0
    X_max = 100
    alpha_min = 0.9
    alpha_max = 1.0
    pa_min = 0.05
    pa_max = 0.25
    N_nest = 25
    N_itertotal = 100
    gamma = 0.1          # noise factor
    lamda = 1.5          # constant used in le'vy flight
    
    total_nodes = 0
    anchor_percent = 0
    t_range = 0 # transmission range
    
    X_anchor = []
    X_unknown = []
    
    X_loc = []

    
    def __init__(self, total_nodes=100, anchor_percent = 0.20, t_range = 25):
        self.total_nodes = total_nodes
        self.anchor_percent = anchor_percent
        self.t_range = t_range
        self.M = int(self.anchor_percent * self.total_nodes) # no. of anchor nodes
        self.N = self.total_nodes - self.M # no. of unknown nodes
        
        
        for i in range(self.M):
            x_anchor = np.random.randint(100)
            y_anchor = np.random.randint(100)
            self.X_anchor.append([x_anchor, y_anchor])
        
        for i in range(self.N):
            x_unknown = np.random.randint(100)
            y_unknown = np.random.randint(100)
            self.X_unknown.append([x_unknown, y_unknown])
        
        self.X_unknown = np.array(self.X_unknown)
        self.X_anchor = np.array(self.X_anchor)
    
        self.X_anchor_og = self.X_anchor.copy()
        self.X_unknown_og = self.X_unknown.copy()
    
    def show_field(self):
        plt.figure(figsize=(8, 8))
        plt.plot(self.X_unknown[:, 0], self.X_unknown[:, 1], 'ro', label="Unknown Node")
        plt.plot(self.X_anchor[:, 0], self.X_anchor[:, 1], 'go', label="Anchor Node")
        plt.legend()
        plt.grid()
        plt.show()
        
        
    # step size (alpha)
    def alpha(self, n_iter):
        return self.alpha_max - ((n_iter/ self.N_itertotal) * (self.alpha_max - self.alpha_min)) # returns step size value

    # Le'vy flight function
    def levy(self):
        temp = np.power(((gamma_function(1 + self.lamda) * np.sin(np.pi * (self.lamda /2))) / (gamma_function((1 + self.lamda)/2) * self.lamda * np.power(2, ((self.lamda - 1)/2)) )), 1/self.lamda)
        u = np.random.normal(0, temp)
        v = np.random.normal(0,1)
        r = u / (np.power(abs(v), (1/self.lamda)))

        return r  # random walk value

    
    # location limit tester
    def limiter(self, point):
        x = point[0]
        y = point[1]
        if x > self.X_max and y > self.X_max:
            x,y = self.X_max, self.X_max
            # X_j = X_rand
        elif x > self.X_max and self.X_min < y < self.X_max:
            x,y = self.X_max, y
            # X_j = X_rand
        elif x > self.X_max and y < self.X_min:
            x,y = self.X_max, self.X_min
            # X_j = X_rand
        elif self.X_min < x < self.X_max and y < self.X_min:
            x,y = x, self.X_min
            # X_j = X_rand
        elif x < self.X_min and y < self.X_min:
            x,y = self.X_min, self.X_min
            # X_j = X_rand
        elif x < self.X_min and self.X_min < y < self.X_max:
            x,y = self.X_min, y
            # X_j = X_rand
        elif x < self.X_min and y > self.X_max:
            x,y = self.X_min, self.X_max
            # X_j = X_rand
        elif self.X_min < x < self.X_max and y > self.X_max:
            x,y = x, self.X_max
            # X_j = X_rand

        return [x,y]

    
    def neighbours(self, node, anchors):
        x = node[0]
        y = node[1]
        X_anchor = anchors
                
        l = []
        for j in range(len(X_anchor)): # for every anchor nodes
            dist_real = np.power((np.power((x - X_anchor[j][0]), 2) + np.power((y - X_anchor[j][1]), 2)), 0.5)
            dist_err = dist_real + np.random.normal(0, (self.gamma*dist_real))

            if dist_err < self.t_range:
                l.append(X_anchor[j])
                
        return l # neighbouring anchors coordinates

    def objective(self, node, n_anchors):
        """objective function (to minimize)"""
        x = node[0]
        y = node[1]
        
        l = self.neighbours(node, n_anchors)
        m = len(l)
        rerror = []
        if len(l) >= 3:
            for ancn in l:
                dist_real = np.power((np.power((x - ancn[0]), 2) + np.power((y - ancn[1]), 2)), 0.5)
                dist_err = dist_real + np.random.normal(0, (self.gamma*dist_real))
                rerror.append(np.power(dist_real - dist_err,2))

        ans = None
        if math.isnan(np.sum(rerror)/m): 
            ans = np.inf
        else:
            ans = np.sum(rerror)/m or None
            
        return ans # mean of square of ranging error

    
    def mod_cs(self, N_anchor):
        X_nest = []
        for i in range(self.N_nest):
            x_nest = np.random.randint(100)
            y_nest = np.random.randint(100)
            X_nest.append([x_nest, y_nest])
        
        Obj_X_nest = []
        for i in range(len(X_nest)):
            Obj_X_nest.append(self.objective(X_nest[i], N_anchor))
        
#         print(Obj_X_nest)
        
        N_iter = 0
        fmins = []
        while(N_iter < self.N_itertotal):
            N_iter += 1
            X_js = []
            for i in range(len(X_nest)):
                X_j = X_nest[i].copy()
                X_j[0] = X_j[0] + self.alpha(N_iter) * self.levy()
                X_j[1] = X_j[1] + self.alpha(N_iter) * self.levy()
                
                X_j = self.limiter(X_j)
                
                F_j = self.objective(X_j, N_anchor) or np.inf
                
                rand_k = np.random.randint(0, len(X_nest))
                
                F_k = self.objective(X_nest[rand_k], N_anchor) or np.inf
                
                if F_j > F_k:
                    X_j[0] = X_nest[rand_k][0]
                    X_j[1] = X_nest[rand_k][1]
                    F_j = F_k

                X_js.append(X_j)
                
            Obj_X_js = []
            for i in range(len(X_js)):
                Obj_X_js.append(self.objective(X_js[i], N_anchor))
            
            Obj_X_js = np.array([np.inf if i is None else i for i in Obj_X_js])
            F_min = Obj_X_js[np.argmin(Obj_X_js)]
#             print(F_min)
            best_sol = X_js[np.argmin(Obj_X_js)]
            
            Pa_j = []
            
            for i in Obj_X_js:
                K = i - F_min
                if K < 1:
                    Pa_j.append(self.pa_min + (self.pa_max - self.pa_min) * K)
                else:
                    Pa_j.append(self.pa_max / N_iter)

            
            for i in range(len(Pa_j)):
                rand_temp = np.random.uniform(0, 1)
                if rand_temp < Pa_j[i]:
                    my_x = np.random.randint(100)
                    my_y = np.random.randint(100)
                    if (self.objective(X_js[i], N_anchor) or np.inf) > (self.objective([my_x, my_y], N_anchor) or np.inf):
                        X_js[i] = [my_x, my_y]
            
            
            Obj_X_js = []
            for i in range(len(X_js)):
                Obj_X_js.append(self.objective(X_js[i], N_anchor))
            
            Obj_X_js = np.array([np.inf if i is None else i for i in Obj_X_js])
            F_min = Obj_X_js[np.argmin(Obj_X_js)]
            print(F_min)
            fmins.append(F_min)
            best_sol = X_js[np.argmin(Obj_X_js)]
            X_nest = X_js.copy()
        
        return best_sol
        
        plt.plot(fmins)
        plt.show()

                
                
    
    def update_Unknown(self, indexes):        
        """Updating Unknown List"""
        
        X_unknown_temp = []
        for j in range(len(self.X_unknown)):
            if j in indexes:
                pass
            else:
                X_unknown_temp.append(self.X_unknown[j])

        self.X_unknown = np.array(X_unknown_temp)
        
    
    
    
    def main(self):
        
        localised_indexes = []
        for i in range(len(self.X_unknown)):
            nn = self.neighbours(self.X_unknown[i], self.X_anchor)
            print(i, end="\r")
            if len(nn) >= 3:
                updated_node_location = self.mod_cs(nn)
                
                self.X_loc.append([updated_node_location, self.X_unknown[i]])
                
                # Updated 
                X_anchor_temp = list(self.X_anchor)
                X_anchor_temp.append(updated_node_location)
                self.X_anchor = np.array(X_anchor_temp)
                
                localised_indexes.append(i)
        
        self.update_Unknown(localised_indexes)
                

                
            
            
        


# In[ ]:





# In[3]:


coco = Cuckoo(anchor_percent=0.1, t_range=25)


# In[ ]:





# In[4]:


coco.main()


# In[5]:


coco.X_loc


# In[7]:


distances = []
for loc in coco.X_loc:
    modfied = loc[0]
    original = loc[1]
    
    distances.append(np.sqrt((modfied[0] - original[0])**2 + (modfied[1] - original[1])**2))


# In[8]:


np.sum(distances)/len(distances)


# In[10]:


len(coco.F_min)


# In[ ]:





# In[129]:


x = [0.5586032962287517,0.6017524859504793,0.7813999837819767,0.46606363790524846,1.0247709434914343,1.252005912323164,0.7495009374755108,0.561087541248166,0.6173482423311187,0.7023584934592145,0.5555991502481207,0.2347329658872088,0.7689648691600328,0.42170816132227973,0.25106692688462434,1.015593603834578,0.409359463643792,0.23540883134280907,0.5913839305269655,0.600203577388256,0.5000388726258884,0.47597790189065975,0.3658604050896738,0.49870120762058423,1.0246665774750374,0.22224505802313932,0.48748272031646006,0.2721460419888341,0.38653619800686356,0.2615371220932912,0.6958970282936564,0.36899753381356065,0.2039305516189484,0.17151067125529823,0.48919234541387874,0.39500189054510565,0.3530413033536565,0.5591781993082433,0.7749257166059389,0.3368538823421103,0.11406348593214892,0.1550569530750901,0.1974787800028229,0.3971344522252228,0.5281760365228262,0.27339861412968086,0.21166372780804943,0.5555694651180548,0.24039256191025846,0.2507917463528349,0.26165577892959974,0.2902138425532703,0.09274393957755404,0.14342988579925614,0.09482794773206088,0.10777099278167086,0.2167607678878698,0.08407868156026704,0.4162032579761468,0.19392965695433126,0.1531513524666788,0.10250629296185206,0.10944110289909664,0.3372574640392889,0.09181169924020706,0.16634369946602098,0.461380034902016,0.10096995500968431,0.5163841659356835,0.4100570120780608,0.3305394441703243,0.12650693376928016,0.4605351112395287,0.4704711004480148,0.39490322471822437,0.18080358326668763,0.40822830771771684,0.23356350554663777,0.22729262817308193,0.24627316532161503,0.16842430253528257,0.4886667259598101,0.17134357949552495,0.3237412112321077,0.21311484660525568,0.40680583941819976,0.2780304862894268,0.31027772973891987,0.07095677497515178,0.2835606346307872,0.19350661176946501,0.2259702437046537,0.23069606587795344,0.10291948260618343,0.3505927932052159,0.520609285578258,0.552334894422106,0.34428465409389686,0.20672654684411254,0.1638661522394548,]


# In[131]:


plt.plot(x)
plt.show()


# In[160]:


coco.X_anchor


# In[161]:





# In[ ]:




