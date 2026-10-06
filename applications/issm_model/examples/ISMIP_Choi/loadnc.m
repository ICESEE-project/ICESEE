function [s]=loadnc(ncfile);
%    loadnc -- function to read the variables from an nc-format file
%
%    usage:
%        s=loadnc(ncfile);
%    where:
%        ncfile     name of nc-format file
%

	if ~exist('ncfile','var') || isempty(ncfile)
		help loadnc
		error('loadnc usage error.');
	end

	[pathstr,name,ext]=fileparts(ncfile);
	if isempty(ext)
		ext='.nc';
	end
	ncfile=fullfile(pathstr,[name ext]);

	a=ncinfo(ncfile);
	disp(sprintf('nc-format file ''%s'' read.',ncfile));

	for i=1:length(a.Variables)
		% matlab reads the dimensions reversed and matrices transposed from netcdf, so compensate for that
		s.(a.Variables(i).Name)=transpose(ncread(ncfile,a.Variables(i).Name));
		disp(sprintf('field ''%s'' of class ''%s'' and size [%dx%d] read.',...
		             a.Variables(i).Name,class(s.(a.Variables(i).Name)),size(s.(a.Variables(i).Name),1),size(s.(a.Variables(i).Name),2)));
	end

end

