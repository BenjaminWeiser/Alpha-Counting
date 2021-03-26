function [Output] = findranges(Act,Ei)
%actinide X is range of alphas, Y is the eloss/depth. Ei is the initial energy
%finds the range at wich the elecctron will com out of the copper with more
%than Eremainung

%E_at_depth give array with depth cordinate and emitted energy of the
%alpha

Eremaining=1.5E6;

X=Act(1:end,1); %Range 
Y=Act(1:end,2); %eloss/depth


Elost=0;
Output= zeros(50, 3);
for i=1:(length(X)-1)
    trap=(Y(i)+Y(i+1))./(2)*(X(i+1)-X(i));
    Elost=Elost+trap;
    Er=Ei-Elost; %remaining E at x
    Output(i,1)=X(i); %depth in column1
    Output(i,2)= Er; %alpha emmited E at that depth in column 2
    if Er<Eremaining
        Output(1,3)=X(i); %range in column 3
    break
    end
end
Output=Output(1:(i-1),:); %trims output to have only nonzero values for depth and E and takes one more off as we want oonly greater than 1.5MeV

end

