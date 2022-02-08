%program to find the effective indexes of guided modes
%in an asymmetric step waveguide

%clear all;
hold on

k=2*pi/.633e-6; 
n1=1.0; 
n2=1.5095; %n2=1.5095; 
n3=1.4711; %n3=1.4711;


nup=n2-0.00001; nlo=n3+0.00001;
for mode_num=0:1:7
for d=0.0:0.01:10
    d=d*1e-6;
    nup=n2-0.00001; nlo=n3+0.00001; %the upper and lower limits to search
    split=(nup+nlo)/2;
    while abs(split-nup)>0.00001
    temp=asym_func_neff(nup,d,k,n1,n2,n3,mode_num)*asym_func_neff(split,d,k,n1,n2,n3,mode_num);
    if temp > 0
        nup=split; 
    else
        nlo=split;
    end
    split=(nup+nlo)/2;
    end
    plot(d*1e6,split,'.r'); hold on
    ylim([n3+0.0001 n2])
end
end
hold off
xlabel('guide thickness (\mum)','Fontsize',20)
ylabel('effective index','Fontsize',20)
gr=gca;
gr.FontSize=20;
    
    
