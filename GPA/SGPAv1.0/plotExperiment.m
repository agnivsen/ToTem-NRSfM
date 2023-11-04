% PLOTEXPERIMENT --- Plot experiment results obtained with genExperiment
% function.
% 
% plotExperiment(name,fignumber,figname,markervector,ColorPalete,leg)
%
% INPUTS
%   
%  name => name of the file that contains a experiment.
%  fignumber => number of figure to display
%  figname => name of the output figure 
%  markervector => cell array of markers used in the graph. by default the
%                   following cell is used {'square','^','diamond','o','p','>','x'};
%
%  ColorPalete => nx3 array with color palete used in the graph. Each row
%                   is a color in [r,g,b] 0<=r,g,b<=1
%                 Default ColorPalete=hsv(nmethods). where nmethods is the
%                 number of methods tested in the selected experiment
% OUTPUTS
% 
% This function has no console outputs. 

% Copyright (C)2010  Adrien Bartoli (1), Daniel Pizarro (2) and Marco Loog
% (3)
% (1) Clermont Université - Clermont-Ferrand, France.
% (2) Universidad de Alcala- Alcala de Henares, Spain
% (3) Delft University - Nederlands
% Send bug reports and feedback to: dani.pizarro@gmail.com
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function plotExperiment(name,fignumber,figname,markervector,ColorPalete,leg)

if(nargin<6)
    if(nargin<5)
    if(nargin<4)
    if(nargin<3)
        if(nargin<2)
            fignumber=1;
            if(nargin<1)
                name=[];
            end
            
        end
        figname='default';
    end
    
        markervector={'square','^','diamond','o','p','>','x'};
  end
ColorPalete=hsv(9);
end
leg=1;
end

  linestyle={'-','--','-.',':','--','-'};
  fontsize=30;
   if(length(name)==0)
  [name,PathName,FilterIndex] = uigetfile('mats/*.mat');
  end
  eval(sprintf('load mats/%s',name));
  timevsparam=zeros(nmethods,nvalues*N);
  errorvsparam=zeros(nmethods,nvalues*N);
nvalues=nvalues;
  if(length(results{nvalues})==0)
nvalues=nvalues-1;
end

  index=1;
  for valueindex=1:nvalues
      for k=[1:N]

          for j=1:nmethods
              resultpervalue=results{valueindex};
              T=resultpervalue{k,j};
              errorvsparam(j,index)=sqrt(T.error);
              timevsparam(j,index)=T.time(1);
          end
          index=index+1;
      end
  end

  hfig=figure(fignumber);

  %ColorPalete=hsv(nmethods);
  for j=[1:nmethods]
      error=zeros(N,nvalues);
      error(:)=errorvsparam(j,1:nvalues*N);
      h=prettyplot(options.parvalues(1:nvalues),error,ColorPalete(j,:),markervector{rem(j,length(markervector))},linestyle{rem(j,length(linestyle))});
      hold on;
      disp(sprintf('Method: %s Min Time: %f  Avgn Time: %f, Max Time %f \n',options.methods{j},min(timevsparam(j,1:nvalues*N)),mean(timevsparam(j,:)),max(timevsparam(j,:))));
      %fprintf('%s & %f |%f |%f \n',options.methods{j},min(timevsparam(j,1:nvalues*N)),mean(timevsparam(j,:)),max(timevsparam(j,:)))  
  end
  h2=get(hfig,'CurrentAxes');
  set(h2,'FontSize',round(fontsize/1.3));
  hold off;
  if(leg==1)
  h=legend(options.methods);
  set(h,'FontSize',fontsize,'Location','Best');
  end
  
  title(sprintf('RMS vs %s',options.varpar),'FontSize',fontsize*1.5);
  
saveas(hfig,figname,'fig');
print(hfig,'-depsc',figname);



