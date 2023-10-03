% Batch-recreate plots 'p--'/'P--' from data saved in 'Rec--' m-files.
% Also updates strings saved in 'Rec--'.
%
%	Loads Rec-files and runs Plots.m on them.
%	Use if you want to change plots and/or name strings without touching
%	underlying data.
%
%   Copy this script to your working folder.
%	Make sure that Plots.m is in Matlab path.
%	Old files are safely moved to _backup_

clear

search = '';
old = {};
new = {}; % (processed sequentially)


flist = dir([pwd '/Rec*' search '*'])
mkdir('_backup_');

for i=1:length(flist)
	
	fprintf('\n%d/%d\n',i,length(flist))
	disp(flist(i).name)
	
	clearvars -except old new flist i
	load(flist(i).name)
	
	for s=1:length(old)
		str1 = strrep(str1, old{s}, new{s});
 		str2 = strrep(str2, old{s}, new{s});
 		str3 = strrep(str3, old{s}, new{s});
	end
	
	% remove old files first, then recalculate
	movefile(flist(i).name, '_backup_');
	flist2 = dir([pwd '/*' flist(i).name(4:end-4) '.png']);
	for j=1:length(flist2)
		movefile(flist2(j).name, '_backup_');
	end 
	
	Plots	
	
end
