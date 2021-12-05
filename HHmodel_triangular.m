close all;clear all;clc;
%============================================================================================================
% Parameters
Cm = 1;%Membrane capacitance (μF/cm^2)
E_Na = 115;%Equilibrium voltage sodium (mV)
E_K=-12; %Equilibrium voltage pottasium (mV)
E_Cl=10.6; %Equilibrium voltage chloride (mV)
G_Na = 120;%Normalization constant Maximum possible conductance mS/cm^2
G_K=36;%mS/cm^2
G_Cl=0.3; % mS/cm^2
%==========================================================================================================
%                                          Determining the Initial values
%===========================================================================================================
Vrest=-65;%Voltage at rest(mV)
n(1) = 0.3;
m(1)=0.04;
h(1)=0.5961;
%=========================================================================================================
%                                  Defining the time vector for simulation
%===========================================================================================================
dt = 0.01; %ms
Tf = 25; %ms
T= 0:dt:Tf;
l=length(T);
%============================================================================================================
%                                               Defining the Current
%===============================================================================================================
I = [0:1:10]; % uA/cm^2
pw = 1; %ms
%===========================================================================================================
%                                Computing coefficients, currents, and derivates 
%===========================================================================================================
for  k = 1 : length(I)
    Vm(k,1) = Vrest;
    for i=1:l
    
        if T(i)>0 && T(i)<=pw
            Iext(k,i)=T(i)*I(k);
        else
            Iext(k,i)=0;
        end
   
    %---calculate the coefficients---%

    %Equations constants
        alpha_n=(0.01*(-55-Vm(k,i)))/(exp( (-55-Vm(k,i)) /10) -1);
        beta_n=0.125*exp( (-65-Vm(k,i))/80 );
        
        alpha_m=(0.1*(-40-Vm(k,i)))/(exp( (-40-Vm(k,i)) /10)-1);
        beta_m=4*exp( (-65-Vm(k,i)) /18);
        
        alpha_h=0.07*exp( (-65-Vm(k,i)) /20); 
        beta_h=1./(exp( (-35-Vm(k,i)) /10)+1);

    %---calculate the derivatives using Euler first order approximation---%
        n(i+1) = n(i) + dt*(alpha_n*(1-n(i)) - beta_n* n(i));
        m(i+1) = m(i) + dt*(alpha_m*(1-m(i)) - beta_m* m(i)); 
        h(i+1) = h(i) + dt*(alpha_h*(1-h(i)) - beta_h* h(i));

    %---calculate Ion conductances---%
        GNa= G_Na*m(i).^3.*h(i);
        gNa(i) = GNa;
        GK= G_K*n(i).^4;
        gK(i) = GK;
     %---calculate the currents---%
        I_Na(i)= GNa*(Vm(k,i)-E_Na-Vrest); 
        I_K(i)= GK*(Vm(k,i)-E_K-Vrest);
        I_Cl(i)= G_Cl*(Vm(k,i)-Vrest-E_Cl);
        Iion = I_K(i) + I_Na(i) +I_Cl(i);
        IION(i) = Iion;
        if T(i)~=Tf   
        Vm(k,i+1)= Vm(k,i) + dt * ((1/Cm)*(Iext(k,i)-Iion));
        end
    end
end
%===========================================================================================================
%                                          Graphs
%===========================================================================================================
figure;
hold on
plot(T,Vm,"LineWidth",2)
legend({'Voltage'})
ylabel('Membrane Voltage (mV)')
xlabel('Time (ms)')
title('Volatge over time in simulated neuron');
hold off

figure;
hold on
plot(T, Iext,"LineWidth",2)
legend({'Extracellular Current'})
ylabel('Pulse Amplitude')
xlabel('Pulse duration')
title('Triangular Monophasic Pulse for stimulation');
hold off

figure;
hold on
plot(T, I_Cl,T, I_Na,T, I_K,T, IION,"LineWidth",2)
legend({'Cl','Na','K','Ionic'})
ylabel('Ionic Current Density (uA/cm)')
xlabel('Time (ms)')
title('Currents');
hold off

figure;
hold on
plot(T,gNa,T,gK,"LineWidth",2)
legend({'Na','K'})
ylabel('Ionic Conductance (ms/cm)')
xlabel('Time (ms)')
title('Conductances');
hold off