%ADI-FDTD for in 3d computational domain with mLor dispersive model

%PEC is @ the boundaries

%Cavity filled with plasma

% K. Prokopidis - 11.5.2017

clear
m0=pi*4e-7; %permability of free space
e0=8.854187817620E-12; %permittivity of free space
c0=1/sqrt(e0*m0); %velocity of light

%dimensions of the cavity
Lx=0.05;
Ly=0.04;
Lz=0.03;

%dimensions of the computational domain
NX=50;
NY=40;
NZ=30; 

CFLN=8; %CFL number
Nt=ceil(5000/CFLN); %time steps

%space step
dx=Lx/(NX-1);
dy=Ly/(NY-1);
dz=Lz/(NZ-1); 

dt=CFLN/(c0*sqrt(1/dx^2+1/dz^2+1/dz^2));  %time step

% Parameters of Drude model (plasma)
omegap=7e9;
vu=0.1*omegap;
epsilon_inf=1;

% Parameters of the mLor model
omega_D=omegap;
gamma=vu;
  
%FDTD coefficients
a_0=e0*omegap^2;
a_1=0;
b_0=0;
b_1=vu;
b_2=1;

%ADI FDTD parameters
%of the dispersive medium
f_1=2./(b_1+(4/dt)*b_2+(b_0/4)*dt);
f_2=-(b_1-(4/dt)*b_2+(b_0/4)*dt)./(b_1+(4/dt)*b_2+(b_0/4)*dt);
f_3=(a_1+(a_0/4)*dt)./(b_1+(4/dt)*b_2+(b_0/4)*dt);
f_4=(-a_1+(a_0/4)*dt)./(b_1+(4/dt)*b_2+(b_0/4)*dt);

h_1=b_1-(4/dt).*b_2;
h_2=b_1+(4/dt).*b_2;

%ADI FDTD parameters for the dispersive medium

%parameters of the electric components
c1x=-dt^2/(4*m0*dx^2);
c1y=-dt^2/(4*m0*dy^2);
c1z=-dt^2/(4*m0*dz^2);

c2x=e0*epsilon_inf+dt^2/(2*m0*dx^2)+sum(f_3);
c2y=e0*epsilon_inf+dt^2/(2*m0*dy^2)+sum(f_3);
c2z=e0*epsilon_inf+dt^2/(2*m0*dz^2)+sum(f_3);

c3=e0*epsilon_inf-sum(f_4);

c4x=dt/(2*dx);
c4y=dt/(2*dy);
c4z=dt/(2*dz);

c5yz=dt^2/(4*m0*dy*dz);
c5xy=dt^2/(4*m0*dx*dy);
c5zx=dt^2/(4*m0*dz*dx);

%parameters of the magnetic components
cHx1=dt/(2*m0*dz);
cHx2=dt/(2*m0*dy);

cHy1=dt/(2*m0*dx);
cHy2=dt/(2*m0*dz);

cHz1=dt/(2*m0*dy);
cHz2=dt/(2*m0*dx);

%source (gaussian)
f_0=8e9;
tw=2/(pi*f_0);
t0=4*tw;
func(1:Nt)=exp(-((1:Nt+0.5)*dt-t0).^2./(tw)^2); 

%Source position
is=NX/2;
js=NY/2;
ks=NZ/2;

%electromagetic fields Ex,Ey,Ez and Hx,Hy,Hz
EX=zeros(NX,NY,NZ);
EX_pre=zeros(NX,NY,NZ); %previous time step

EY=zeros(NX,NY,NZ);
EY_pre=zeros(NX,NY,NZ); %previous time step

EZ=zeros(NX,NY,NZ);
EZ_pre=zeros(NX,NY,NZ); %previous time step

HX=zeros(NX-1,NY-1,NZ-1);
HY=zeros(NX-1,NY-1,NZ-1);
HZ=zeros(NX-1,NY-1,NZ-1);

%auxiliary variables
PX=zeros(NX,NY,NZ);
temp_PX=zeros(NX,NY,NZ);
QX=zeros(NX,NY,NZ);

PY=zeros(NX,NY,NZ);
temp_PY=zeros(NX,NY,NZ);
QY=zeros(NX,NY,NZ);

PZ=zeros(NX,NY,NZ);
temp_PZ=zeros(NX,NY,NZ);
QZ=zeros(NX,NY,NZ);

%--------------- Build the matrix A ---------------
A_z=sparse(NZ,NZ); %for the update of EX at first substep abd EY at second substep
A_y=sparse(NY,NY); %for the update of EZ at first substep and EX at second substep
A_x=sparse(NX,NX); %for the update of EY at first substup and EZ at second substep

%matrix b_x, b_y and b_z
b_x=zeros(NX,1);
b_y=zeros(NY,1);
b_z=zeros(NZ,1);

%Observation point
EZ_ADI=zeros(Nt,1);

%First we build the sparse arrays A_x,A_y, and A_z which are constant during the
%simulation - we assume homogeneous medium

A_z(1,1)=1; %PEC
for k=2:NZ-1        
                 
        %array A_z
        A_z(k,k-1)=c1z;
        A_z(k,k)=c2z;
        A_z(k,k+1)=c1z;
end
A_z(NZ,NZ)=1;  %PEC

b_z(1)=0; %PEC
b_z(NZ)=0; %PEC


A_x(1,1)=1; %PEC
for i=2:NX-1        
                 
        %array A_x
        A_x(i,i-1)=c1x;
        A_x(i,i)=c2x;
        A_x(i,i+1)=c1x;
end
A_x(NX,NX)=1;  %PEC

b_x(1)=0; %PEC
b_x(NX)=0; %PEC

A_y(1,1)=1; %PEC
for j=2:NY-1        
                 
        %array A_y
        A_y(j,j-1)=c1y;
        A_y(j,j)=c2y;
        A_y(j,j+1)=c1y;
end
A_y(NY,NY)=1;  %PEC

b_y(1)=0; %PEC
b_y(NY)=0; %PEC

A_x_inv=inv(A_x);
A_y_inv=inv(A_y);
A_z_inv=inv(A_z);


tic;
%----------------- ADI-FDTD simulation ----------------------------------
for n=1:Nt
    disp(n);
    
    % -------------------------- First substep (n -> n+1/2) ----------------------------------
       
    
   % ------------------- EX update ----------------- (implicitly)
   %saves past times of EX @ n
    EX_pre=EX(1:NX,1:NY,1:NZ);
   
   % a linear system of equations is solved for every i and j
   for i=1:NX-1
       for j=2:NY-1
              
       %A_z(1,1)=1; %PEC
       %A_z(NZ,NZ)=1;  %PEC
       
       % --------------- build matrix b
       
       for k=2:NZ-1        
                 
        %array A_z
        %A_z(k,k-1)=c1z;
        %A_z(k,k)=c2z;
        %A_z(k,k+1)=c1z;
            
        %column b            
        
         b_z(k)=c3*EX(i,j,k)-f_1*QX(i,j,k)+(1-f_2)*PX(i,j,k)...
             +c4y*(HZ(i,j,k)-HZ(i,j-1,k))...
             -c4z*(HY(i,j,k)-HY(i,j,k-1))...
             -c5zx*(EZ(i+1,j,k)-EZ(i,j,k)-EZ(i+1,j,k-1)+EZ(i,j,k-1));
        end       
         
       %find EX implicitly solving the linear system
       %solution using Thomas algorithm
       %EX(i,j,1:NZ)=A_z\b_z;
       EX(i,j,1:NZ)=A_z_inv*b_z;
       end
    
   end
   
    % ---------------- EY update -------------------- (implicitly)
     %saves past times of EY @ n
    EY_pre=EY(1:NX,1:NY,1:NZ);
    
    % a linear system of equations is solved for every j and k
    for j=1:NY-1
        for k=2:NZ-1
              
       %A_x(1,1)=1; %PEC
       %A_x(NX,NX)=1;  %PEC
       
       % --------------- build matrix b
       
       for i=2:NX-1        
                 
%         %array A_x
%         A_x(i,i-1)=c1x;
%         A_x(i,i)=c2x;
%         A_x(i,i+1)=c1x;
            
        %column b            
                            
         b_x(i)=c3*EY(i,j,k)-f_1*QY(i,j,k)+(1-f_2)*PY(i,j,k)...
             +c4z*(HX(i,j,k)-HX(i,j,k-1))...
             -c4x*(HZ(i,j,k)-HZ(i-1,j,k))...
             -c5xy*(EX_pre(i,j+1,k)-EX_pre(i,j,k)-EX_pre(i-1,j+1,k)+EX_pre(i-1,j,k)); 
       end   
       
    
       %find EY implicitly solving the linear system
       %solution using Thomas algorithm
       %EY(1:NX,j,k)=A_x\b_x;
       EY(1:NX,j,k)=A_x_inv*b_x;
       end
    
   end 
    
    
    % --------------- EZ update --------------------- (implicitly)
    %saves past times of EZ @ n
    EZ_pre=EZ(1:NX,1:NY,1:NZ);
    
     % a linear system of equations is solved for every i and j
   for i=2:NX-1
       for k=1:NZ-1
              
       %A_y(1,1)=1; %PEC
       %A_y(NY,NY)=1;  %PEC
       
       % --------------- build matrix b
       
       for j=2:NY-1        
                 
        %array A_y
        %A_y(j,j-1)=c1y;
        %A_y(j,j)=c2y;
        %A_y(j,j+1)=c1y;
            
        %column b            
        
         if (i==is && j==js && k==ks)
             
         b_y(j)=c3*EZ(i,j,k)-f_1*QZ(i,j,k)+(1-f_2)*PZ(i,j,k)...
             +c4x*(HY(i,j,k)-HY(i-1,j,k))...
             -c4y*(HX(i,j,k)-HX(i,j-1,k))...
             -c5yz*(EY_pre(i,j,k+1)-EY_pre(i,j,k)-EY_pre(i,j-1,k+1)+EY_pre(i,j-1,k))-(0.5*dt/dx/dy)*func(n); 
             
             
         else 
         b_y(j)=c3*EZ(i,j,k)-f_1*QZ(i,j,k)+(1-f_2)*PZ(i,j,k)...
             +c4x*(HY(i,j,k)-HY(i-1,j,k))...
             -c4y*(HX(i,j,k)-HX(i,j-1,k))...
             -c5yz*(EY_pre(i,j,k+1)-EY_pre(i,j,k)-EY_pre(i,j-1,k+1)+EY_pre(i,j-1,k)); 
        end
       end    
       
    
        %find EZ implicitly solving the linear system
        %solution using Thomas algorithm
        %EZ(i,1:NY,k)=A_y\b_y;
        EZ(i,1:NY,k)=A_y_inv*b_y;
        end
    
   end
   % -----------------  source excitation
   %EZ(is,js,ks)=EZ(is,js,ks)-0.5*(dt/e0/epsilon_inf/dx/dy)*func(n);   
   
   
    % --------------- Update magnetic fields (explicitly) 
    %update of HX @ (n+1/2) - explicit 
    HX(1:NX-1,1:NY-1,1:NZ-1)=HX(1:NX-1,1:NY-1,1:NZ-1)...
        +cHx1.*(EY_pre(1:NX-1,1:NY-1,2:NZ)-EY_pre(1:NX-1,1:NY-1,1:NZ-1))...
        -cHx2.*(EZ(1:NX-1,2:NY,1:NZ-1)-EZ(1:NX-1,1:NY-1,1:NZ-1));
    
    
    %update of HY @ (n+1/2) - explicit 
    HY(1:NX-1,1:NY-1,1:NZ-1)=HY(1:NX-1,1:NY-1,1:NZ-1)...
        +cHy1.*(EZ_pre(2:NX,1:NY-1,1:NZ-1)-EZ_pre(1:NX-1,1:NY-1,1:NZ-1))...
        -cHy2.*(EX(1:NX-1,1:NY-1,2:NZ)-EX(1:NX-1,1:NY-1,1:NZ-1));
    
    %update of HZ @ (n+1/2) - explicit
    HZ(1:NX-1,1:NY-1,1:NZ-1)=HZ(1:NX-1,1:NY-1,1:NZ-1)...
        +cHz1.*(EX_pre(1:NX-1,2:NY,1:NZ-1)-EX_pre(1:NX-1,1:NY-1,1:NZ-1))...
        -cHz2.*(EY(2:NX,1:NY-1,1:NZ-1)-EY(1:NX-1,1:NY-1,1:NZ-1));
    
    
    % --------------- Update auxiliary variables P and Q                               
    
    %update Px
    
       temp_PX(1:NX,1:NY,1:NZ)=PX(1:NX,1:NY,1:NZ);
        
       PX(1:NX,1:NY,1:NZ)=f_1*QX(1:NX,1:NY,1:NZ)+f_2*temp_PX(1:NX,1:NY,1:NZ)...
           +f_3*EX(1:NX,1:NY,1:NZ)+f_4*EX_pre(1:NX,1:NY,1:NZ);
    
   
    %update Py
    
       temp_PY(1:NX,1:NY,1:NZ)=PY(1:NX,1:NY,1:NZ);
        
       PY(1:NX,1:NY,1:NZ)=f_1*QY(1:NX,1:NY,1:NZ)+f_2*temp_PY(1:NX,1:NY,1:NZ)...
           +f_3*EY(1:NX,1:NY,1:NZ)+f_4*EY_pre(1:NX,1:NY,1:NZ);
    
    
    %update Pz
    
       temp_PZ(1:NX,1:NY,1:NZ)=PZ(1:NX,1:NY,1:NZ);
        
       PZ(1:NX,1:NY,1:NZ)=f_1*QZ(1:NX,1:NY,1:NZ)+f_2*temp_PZ(1:NX,1:NY,1:NZ)...
           +f_3*EZ(1:NX,1:NY,1:NZ)+f_4*EZ_pre(1:NX,1:NY,1:NZ);
    
    
    
    %update Qx    
       QX(1:NX,1:NY,1:NZ)=-QX(1:NX,1:NY,1:NZ)+h_1*temp_PX(1:NX,1:NY,1:NZ)+h_2*PX(1:NX,1:NY,1:NZ);
        
    %update Qy
       QY(1:NX,1:NY,1:NZ)=-QY(1:NX,1:NY,1:NZ)+h_1*temp_PY(1:NX,1:NY,1:NZ)+h_2*PY(1:NX,1:NY,1:NZ);
       
    %update Qz
       QZ(1:NX,1:NY,1:NZ)=-QZ(1:NX,1:NY,1:NZ)+h_1*temp_PZ(1:NX,1:NY,1:NZ)+h_2*PZ(1:NX,1:NY,1:NZ);
    
    
    
    %------------------------- Second substep (n+1/2 -> n) --------------------
 
    
       
    % ------------------- EX update ----------------- (implicitly)
    %save the EX of n+1/2 
    EX_pre=EX(1:NX,1:NY,1:NZ); 
    
    % a linear system of equations is solved for every i and k
   for i=1:NX-1
       for k=2:NZ-1
              
       %A_y(1,1)=1; %PEC
       %A_y(NY,NY)=1;  %PEC
       
       % --------------- build matrix b
       
       for j=2:NY-1        
                 
        %array A_y
        %A_y(j,j-1)=c1y;
        %A_y(j,j)=c2y;
        %A_y(j,j+1)=c1y;
            
        %column b            
        
                     
         b_y(j)=c3*EX(i,j,k)-f_1*QX(i,j,k)+(1-f_2)*PX(i,j,k)...
             +c4y*(HZ(i,j,k)-HZ(i,j-1,k))...
             -c4z*(HY(i,j,k)-HY(i,j,k-1))...
             -c5xy*(EY(i+1,j,k)-EY(i,j,k)-EY(i+1,j-1,k)+EY(i,j-1,k));
       end
        
    
        %find EX implicitly solving the linear system
        %solution using Thomas algorithm
        %EX(i,1:NY,k)=A_y\b_y;
        EX(i,1:NY,k)=A_y_inv*b_y;
        end
    
   end
    
    
    % ------------------- EY update ----------------- (implicitly)
     %saves past times of EY @ n
    EY_pre=EY(1:NX,1:NY,1:NZ);
    
    % a linear system of equations is solved for every j and k
    for i=2:NX-1
        for j=1:NY-1
              
       %A_z(1,1)=1; %PEC
       %A_z(NZ,NZ)=1;  %PEC
       
       % --------------- build matrix b
       
       for k=2:NZ-1        
                 
        %array A_z
        %A_z(k,k-1)=c1z;
        %A_z(k,k)=c2z;
        %A_z(k,k+1)=c1z;
            
        %column b            
        
                     
         b_z(k)=c3*EY(i,j,k)-f_1*QY(i,j,k)+(1-f_2)*PY(i,j,k)...
             +c4z*(HX(i,j,k)-HX(i,j,k-1))...
             -c4x*(HZ(i,j,k)-HZ(i-1,j,k))...
             -c5yz*(EZ(i,j+1,k)-EZ(i,j,k)-EZ(i,j+1,k-1)+EZ(i,j,k-1)); 
       end
       
    
        %find EY implicitly solving the linear system
        %solution using Thomas algorithm
        %EY(i,j,1:NZ)=A_z\b_z;
        EY(i,j,1:NZ)=A_z_inv*b_z;
        end
    
   end 
    
    % ------------------- EZ update ----------------- (implicitly)
    %saves past times of EZ @ n
    EZ_pre=EZ(1:NX,1:NY,1:NZ);
    
     % a linear system of equations is solved for every j and k
   for j=2:NY-1
       for k=1:NZ-1
              
       %A_x(1,1)=1; %PEC
       %A_x(NX,NX)=1;  %PEC
       
       % --------------- build matrix b
       
       for i=2:NX-1        
                 
        %array A_x
        %A_x(i,i-1)=c1x;
        %A_x(i,i)=c2x;
        %A_x(i,i+1)=c1x;
            
        %column b            
        
         
         if (i==is && j==js && k==ks)
          b_x(i)=c3*EZ(i,j,k)-f_1*QZ(i,j,k)+(1-f_2)*PZ(i,j,k)...
             +c4x*(HY(i,j,k)-HY(i-1,j,k))...
             -c4y*(HX(i,j,k)-HX(i,j-1,k))...
             -c5zx*(EX_pre(i,j,k+1)-EX_pre(i,j,k)-EX_pre(i-1,j,k+1)+EX_pre(i-1,j,k))-(0.5*dt/dx/dy)*func(n);   
             
             
        else
         b_x(i)=c3*EZ(i,j,k)-f_1*QZ(i,j,k)+(1-f_2)*PZ(i,j,k)...
             +c4x*(HY(i,j,k)-HY(i-1,j,k))...
             -c4y*(HX(i,j,k)-HX(i,j-1,k))...
             -c5zx*(EX_pre(i,j,k+1)-EX_pre(i,j,k)-EX_pre(i-1,j,k+1)+EX_pre(i-1,j,k));
        end
       end
       
    
        %find EZ implicitly solving the linear system
        %solution using Thomas algorithm
        %EZ(1:NX,j,k)=A_x\b_x;
        EZ(1:NX,j,k)=A_x_inv*b_x;
        end
    
   end  
   % -----------------  source excitation
   %EZ(is,js,ks)=EZ(is,js,ks)-0.5*(dt/e0/epsilon_inf/dx/dy)*func(n);   
    
    % --------------- Update magnetic fields (explicitly) 
     %update of HX @ (n+1) - explicit 
    HX(1:NX-1,1:NY-1,1:NZ-1)=HX(1:NX-1,1:NY-1,1:NZ-1)...
        +cHx1.*(EY(1:NX-1,1:NY-1,2:NZ)-EY(1:NX-1,1:NY-1,1:NZ-1))...
        -cHx2.*(EZ_pre(1:NX-1,2:NY,1:NZ-1)-EZ_pre(1:NX-1,1:NY-1,1:NZ-1));
    
    
    %update of HY @ (n+1) - explicit 
    HY(1:NX-1,1:NY-1,1:NZ-1)=HY(1:NX-1,1:NY-1,1:NZ-1)...
        +cHy1.*(EZ(2:NX,1:NY-1,1:NZ-1)-EZ(1:NX-1,1:NY-1,1:NZ-1))...
        -cHy2.*(EX_pre(1:NX-1,1:NY-1,2:NZ)-EX_pre(1:NX-1,1:NY-1,1:NZ-1));
    
    %update of HZ @ (n+1) - explicit
    HZ(1:NX-1,1:NY-1,1:NZ-1)=HZ(1:NX-1,1:NY-1,1:NZ-1)...
        +cHz1.*(EX(1:NX-1,2:NY,1:NZ-1)-EX(1:NX-1,1:NY-1,1:NZ-1))...
        -cHz2.*(EY_pre(2:NX,1:NY-1,1:NZ-1)-EY_pre(1:NX-1,1:NY-1,1:NZ-1));
     
    
                                   
    % --------------- Update auxiliary variables P and Q  
    %update Px
    
       temp_PX(1:NX,1:NY,1:NZ)=PX(1:NX,1:NY,1:NZ);
        
       PX(1:NX,1:NY,1:NZ)=f_1*QX(1:NX,1:NY,1:NZ)+f_2*temp_PX(1:NX,1:NY,1:NZ)...
           +f_3*EX(1:NX,1:NY,1:NZ)+f_4*EX_pre(1:NX,1:NY,1:NZ);
    
   
    %update Py
    
       temp_PY(1:NX,1:NY,1:NZ)=PY(1:NX,1:NY,1:NZ);
        
       PY(1:NX,1:NY,1:NZ)=f_1*QY(1:NX,1:NY,1:NZ)+f_2*temp_PY(1:NX,1:NY,1:NZ)...
           +f_3*EY(1:NX,1:NY,1:NZ)+f_4*EY_pre(1:NX,1:NY,1:NZ);
    
    
    %update Pz
    
       temp_PZ(1:NX,1:NY,1:NZ)=PZ(1:NX,1:NY,1:NZ);
        
       PZ(1:NX,1:NY,1:NZ)=f_1*QZ(1:NX,1:NY,1:NZ)+f_2*temp_PZ(1:NX,1:NY,1:NZ)...
           +f_3*EZ(1:NX,1:NY,1:NZ)+f_4*EZ_pre(1:NX,1:NY,1:NZ);
    
    
    
    %update Qx
    
       QX(1:NX,1:NY,1:NZ)=-QX(1:NX,1:NY,1:NZ)+h_1*temp_PX(1:NX,1:NY,1:NZ)+h_2*PX(1:NX,1:NY,1:NZ);
    
    
    %update Qy
    
       QY(1:NX,1:NY,1:NZ)=-QY(1:NX,1:NY,1:NZ)+h_1*temp_PY(1:NX,1:NY,1:NZ)+h_2*PY(1:NX,1:NY,1:NZ);
    
    
    %update Qz
    
       QZ(1:NX,1:NY,1:NZ)=-QZ(1:NX,1:NY,1:NZ)+h_1*temp_PZ(1:NX,1:NY,1:NZ)+h_2*PZ(1:NX,1:NY,1:NZ);
    
    
    
        
   %Observation point
    EZ_ADI(n)=EZ(4,4,4);
   
  
    
end

elapsed_time_in_sec=toc;
  disp('Elapsed time in min');
  elapsed_time_in_sec/60


t_ADI=dt.*(1:Nt);
 
%saving to file (error norm)
filename=strcat('P_obs_ADI_',num2str(CFLN),'.mat');
save(filename,'dt','Nt','t_ADI','EZ_ADI');

