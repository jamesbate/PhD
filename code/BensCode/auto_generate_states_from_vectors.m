%build optimal sensor states, S qubit case (one at each site in 1D), Taylor expansion
clear;clc

Hv=[1 0]';
Vv=[0,1]'; %H and v are my qubit logic 0 and 1 states, becasue I'm a photon guy at heart 

S=5; %number of sensors/qubits
O=3; %order of spatial terms (to be 'sensed') to go up to in Taylor expansion.

%if mod(S,2)==0  %even

%%positions of ions in dimensionless units
    %m=40; N=5; nu=0.500E6;
    %posvec=IonPositionsLinearString_fun(nu,m,N);
    %posvec=posvec*1E6;posvec/abs(posvec(2)) 
    %
    %   -2.1201   -1.0000   -0.0000    1.0000    2.1201
    X=[-2.12,-1,0,1,2.12];

    %X=-(S-1)/2:1:(S-1)/2; %poisitions centre on zero. spacing for odd qubit number
    %   -3.5000   -2.5000   -1.5000   -0.5000    0.5000    1.5000    2.5000    3.5000
    %else
    %x=-(S)/2:1:(S)/2; %poisitions centre on zero. spacing for even sensor number
    %-1.7500   -1.2500   -0.7500   -0.2500    0.2500    0.7500    1.2500    1.7500
    %end

%centre of Taylor expansion
    XO=0;


%Taylor expansion
    for or=1:(O)
    for jj=1:length(X);
        B(jj,or)=(X(jj)-XO)^(or-1);
        %if X(jj)==3
        %    B(jj,or)=-1*B(jj,or);
        %end
    end
    Bnorm(:,or)=B(:,or)/norm(B(:,or));
    end
    B
    Bnorm


% %some jiggery pokery here, opposite sign sensors every second site. 
% mu=[1,1,1,1;-1,-1,-1,-1;1,1,1,1;-1,-1,-1,-1]
% te=zeros(4,5); te(:,2:5)=mu.*Bnorm; te(:,1)=1; Bnorm=te
% te2=zeros(4,5); te2(:,2:5)=mu.*B; te2(:,1)=1; B=te2
% S=S+1;

%I've forgotten why this is here and what this does
    kernel=null(Bnorm');


% %Fourier series for our cavity, example. 
% k1=854;
% k2=806;
% 
% X=X.*pi/(1*k1);
% for jj=1:length(X);
%     B(jj,1)=sin(k1*X(jj)+pi/2)^2;
%     B(jj,2)=sin(k2*X(jj)+0.9+pi/2)^2;
% end
% Bnorm(:,1)=B(:,1)/norm(B(:,1));
% Bnorm(:,2)=B(:,2)/norm(B(:,2));
% B
% Bnorm


%now find optimal state vectors
    A=zeros(size(Bnorm));


    for hh=1:O
        tem=Bnorm(:,hh);
        for kk=1:O;
            if kk~=hh;
            tem=tem-(Bnorm(:,kk)'*tem)*Bnorm(:,kk);
            end
        end
        A(:,hh)=tem;
        Anorm(:,hh)=tem/max(abs(tem));
    end

    A

%clean up very small non-zero values at 10-16 level.
    Anorm(abs(Anorm)<=1E-10)=0;
    A(abs(A)<=1E-10)=0;


%attempt at obtaining/quantify some kind of sensitivity
    for gg=1:O;
    Sens(gg)=Anorm(:,gg)'*B(:,gg);
    end
    Sens

%sum(abs(Sens))




%%
%build entangled state vectors, according to the Dür receipe learned in his
%office

    Bp = cell(2,S,or);
    Bm = cell(2,S,or);


    for v=1:2; %these are the two terms in superposition
        if v==2
            Anorm=-1.*Anorm; %probably redundant now that use normalised A vectors?
        end

        for h=1:or;
            for g=1:S;
                de=sign(Anorm(g,h));
                if de==1
                    Bp{v,g,h}=Hv;
                    np(v,g,h)=0;
                elseif de==0
                    Bp{v,g,h}=Hv;
                    np(v,g,h)=0;
                elseif de==-1;
                    Bp{v,g,h}=Vv;
                    np(v,g,h)=1;
                end    
            end
        end
    end

    %Bp{term, qubit number, taylor order}

    psivec=zeros(2^S,or);

    for bb=1:or; %loop over number of Taylor orders
        T1=1;T2=1;
        for aa=1:S %nested loop to make S qubit tensor product
        T1=nkron(T1,Bp{1,aa,bb});
        T2=nkron(T2,Bp{2,aa,bb});
        end
        psivec(:,bb)=T1+T2;
    end

%now work out flip schedule in super Ramsey
    flipmat=zeros(S,or);
    for a=1:S;
        for b=1:or;
            flipmat(a,b)=abs(Anorm(a,b))/2+1/2;
        end
    end

    psivec
    flipmat

    A = char.empty;
    emp=zeros(size(psivec));
    for ll=1:size(psivec,2);
    kk=1;
        for ii=1:length(psivec);
            if psivec(ii,ll)==1       
                A{kk,ll}=dec2bin(ii-1);
                kk=kk+1;
            end
        end
    end


      A      
