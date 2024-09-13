import numpy as np

# function to initialize s,c,g in the domain
# flag = 0 is no surfactantimplementation#
# flag = 1 is with surfactant
# 1-s0 = initial residual saturation 
# c0 = concentration of polymer in injected mixture
# g0 = concentration of surfactant in injected mixture

### Vectorized implementation
def s0c0(para_box, phi, s_0, c_0, g_0, flag):
    s0 = np.zeros((para_box["n"] + 1, para_box["m"] + 1))
    c0 = np.copy(s0)
    D = (phi > 1e-10) + (np.abs(phi) < 1e-10)
    
#     if flag ==0
#         g0=[];
#         s0=(1-s_0)*ones(para.box.n+1,para.box.m+1);
#     elseif flag == 1
#         s0=(1-s_0)*ones(para.box.n+1,para.box.m+1);
#         g0=c0;
#     end
#     

    if flag == 0:
        g0 = []
        s0 = (~D) + D * (1 - s_0)
        c0 = (~D) * c_0
    elif flag == 1:
        s0 = (~D) + D * (1 - s_0)
        c0 = (~D) * c_0
        g0 = (~D) * g_0
        
    return s0, c0, g0
        
### For loop implementation #########

#     s0=zeros(para.box.n+1,para.box.m+1);
#     c0=s0;
#     if flag == 0
#         g0 = [];
#         for i=1:para.box.m+1
#             for j=1:para.box.n+1
#                 if phi(j,i)>10^(-10)||abs(phi(j,i))<10^(-10)
#                     s0(j,i)=1-s_0;
#                     c0(j,i)=0;
#                 else
#                     s0(j,i)=1;
#                     c0(j,i)=c_0;
#                 end
#             end
#         end
#     elseif flag == 1    
#         g0=zeros(para.box.n+1,para.box.m+1);
#         for i=1:para.box.m+1
#             for j=1:para.box.n+1
#                 if phi(j,i)>10^(-10)||abs(phi(j,i))<10^(-10)
#                     s0(j,i)=1-s_0;
#                     c0(j,i)=0;
#                     g0(j,i)=0;
#                 else
#                     s0(j,i)=1;
#                     c0(j,i)=c_0;
#                     g0(j,i)=g_0;
#                 end
#             end
#         end
#     end

### Testing of vector implementation
#     isequal(s0,s01)
#     isequal(c0,c01)
#     pause   