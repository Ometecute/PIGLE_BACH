

tmp=split(pwd,'/');
if ~strcmp(tmp{end},'mexInteractions')
    delete *mexa64 frmBuilder* interp1lin* intrAdsr* setclist* *slxc
else
    disp('do not delete core files')
end




if 0 
%%
    tmp = dir('**/*frmBuilder*');
    for i = 1:length(tmp)
        str = [tmp(i).folder,'/',tmp(i).name];
        delete(str)
    end
    tmp = dir('**/*interp1lin*');
    for i = 1:length(tmp)
        str = [tmp(i).folder,'/',tmp(i).name];
        delete(str)
    end
    tmp = dir('**/*intrAdsr*');
    for i = 1:length(tmp)
        str = [tmp(i).folder,'/',tmp(i).name];
        delete(str)
    end
    tmp = dir('**/*setclist*');
    for i = 1:length(tmp)
        str = [tmp(i).folder,'/',tmp(i).name];
        delete(str)
    end
    tmp = dir('**/*slxc');
    for i = 1:length(tmp)
        str = [tmp(i).folder,'/',tmp(i).name];
        delete(str)
    end
end