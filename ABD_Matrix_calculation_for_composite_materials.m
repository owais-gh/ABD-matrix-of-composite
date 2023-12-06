%Classical Laminate Theory 
%Calculating ABD matrix for single material composite with equal thickness layers


tic;                      %--Code run time
clear all; clc; format long g  %clear any existing assigned values|clear command window|Set output display format to decimals

%--Material data--
t=125e-3;                 %Thickness of layers/This code works only for equal thickness layers composite materials
e11=1.35e5;
e22=1e4;
g12=5e3;
v12=0.27;
n=5;                      %--Number of layers|Make sure the number of layes and number of orientations are equal--
theta=[0 90 90 0];        %--Enter the layer orientation in deg|Works both for symmetric and unsymmetric layer orientations--

%if n!=length(theta)
%    display("Number of layers(n) and Number of layer orientations are not equal(theta)|Displayed vlaues are wrong")

zp=linspace(1,n,n);       %Creating a array of the layer positions

%--Calculating Poison's ratio, v21--
v21=v12*e22/e11;

%--reduced matrix for material--
q11=e11/(1-(v12*v21));
q22=e22/(1-(v12*v21));
q12=(v12*e22)/(1-(v12*v21));
q66=g12;

%--Creating stiffness matrix--
q=[q11 q12 0; q12 q22 0; 0 0 q66]; %delete ";" to display


%--z height--
z=[];
z0=-(n*t)/2;
f=0;
for h=1:n+1;
    zh=z0+(f*t);
    z=[z,zh];
    f=f+1;
end

%--Calculating transformed reduced matrix--
a=0;
for k=1:n;
    
    c=cosd(theta(k));
    s=sind(theta(k));
    
    qt(1,1)=(q11*c^4)+(2*(q12+(2*q66))*(c^2)*(s^2))+(q22*(s^4));
    qt(1,2)=((q11+q22-(4*q66))*((c^2)*(s^2)))+(q12*((c^4)+(s^4)));
    qt(1,3)=((q11-q12-(2*q66))*((c^3)*s))+((q12-q22+(2*q66))*(c*(s^3)));
    qt(2,1)=qt(1,2);
    qt(2,2)=(q11*s^4)+(2*(q12+(2*q66))*((c^2)*(s^2)))+(q22*(c^4));
    qt(2,3)=((q11-q12-(2*q66))*(c*(s^3)))+((q12-q22+(2*q66))*((c^3)*s));
    qt(3,1)=qt(1,3);
    qt(3,2)=qt(2,3);
    qt(3,3)=((q11+q22-(2*(q12+q66)))*((c^2)*(s^2)))+(q66*((c^4)+(s^4)));
    
    a=a+1;
    qta(:,:,a)=qt;  %Creating an empty matrix of nilxnilxnil size matrix, (:,:,"matrix index")
    
end

%Displaying transformed matrix for layers
%{
qta(:,:,2);  %--Tranformed matrix for specific layer, eg: 2=layer position
qta;         %--Tranformed matrix for all the layers
%}


%--Calculating A, B & D matrix-- 
for i=1:3;
    for j=1:3;
        for x=1:n;
            Aa(i,j,x)=qta(i,j,x)*(z(x+1)-z(x)); 
            Bb(i,j,x)=(qta(i,j,x)*((z(x+1)^2)-(z(x)^2)))/2; 
            Dd(i,j,x)=(qta(i,j,x)*((z(x+1)^3)-(z(x)^3)))/3; 
        end
    end
end


%--Display A,B&D matrix--
A=sum(Aa,3);           %delete ";" to display
B=sum(Bb,3);           %delete ";" to display
D=sum(Dd,3);           %delete ";" to display

%--Display ABD matrix as 6x6 matrix [A B;B D]
ABD_matrix=[A B;B D];  %delete ";" to display


%--For representing the stiffness values in each layer, stairs graph|Uncomment if needed--
%{
zpl=[0,zp];
Agp=[0];
Bgp=[0];
Dgp=[0];
for f=1:n
    Ag=Aa(1,1,f);
    Bg=Bb(1,1,f);
    Dg=Dd(1,1,f);
    Agp=[Agp,Ag];
    Bgp=[Bgp,Bg];
    Dgp=[Dgp,Dg];
end
Agp;
Bgp;
Dgp;
stairs(repmat(Agp,1,1),zpl,'DisplayName','Memebrane stiffness - A','LineWidth',2)  
hold on
stairs(repmat(Bgp,1,1),zpl,'DisplayName','Coupling stiffness - B','LineWidth',2)
stairs(repmat(Dgp,1,1),zpl,'DisplayName','Plate stiffness - D','LineWidth',2)
hold off
legend
title('Stiffness value for the layer position')
xlabel('Stiffness value')
ylabel('Layer position')
grid on

%view([90 -90]) %--flips the x&y axis for viewing--xaxis=layer position&yaxis=stiffness value
%}

toc  %--Code run time

%--*************************--%