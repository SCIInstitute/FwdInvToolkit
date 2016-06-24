function LOW=lowpassma(M,win)
%function LOW=lowpassma(M,win)
% lowpass Moving Average filter applied to rows of matrix M; finite impulse
% response; no phase shift
% zero output for frequencies at multiples of  f=win*sampling_interval
% for highpass: use PSI=PSI-lowpassma(M,win)
% 2005-02-17;  %V. Jacquemet/ A. van Oosterom
% 2008-02-12; AvO winb=round(   changed into:   winb=floor(  
LOW=M;
if win <2, return, end
trans=0;
[nsig nt]=size(M);
if nt==1,
    trasp=1;
    M=M';
    [nsig nt]=size(M);
end
winb=floor(win/2);
wine=win-winb;
LEAD=M(:,1)*ones(1,winb);
TRAIL=M(:,nt)*ones(1,wine);
M=[LEAD M TRAIL];
X = cumsum([zeros(nsig,1) M],2);
LOW = X(:,win+1:win+nt)-X(:,1:nt);
LOW = LOW / win;
if trans==1,
    LOW=LOW';
end