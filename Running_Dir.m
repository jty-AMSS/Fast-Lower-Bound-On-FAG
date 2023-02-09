% Running_Dir
% Running the file under 'dir'
cd dir;
X=dir();
Name={};
for i=1:length(X)
    if(length(X(i).name())>3)
        Name{end+1}=X(i).name();
    end
end
cd ..
Table=[];
o_path=pwd;
for i=1:length(Name)
    Temp_name=[o_path ,'\dir\', Name{i}];
    disp(Temp_name)
    tic;
    [lb_l,Q,f,X0,Index]=FSOSBulider(Temp_name,d,sparsity);
    time=toc;
    tline=[lb_l,time,d,size(Index,1)];
    Table=[Table;tline];
end

