function M = design_to_mdd_python(file)
% file = 'test_1.csv'

% mooring configuration
opts = detectImportOptions(file,'VariableNamingRule','preserve');
opts = setvartype(opts,{'serial'}, 'string');
mtab  = readtable(file,opts); % load mooring


% clean up table
idx   = cellfun(@(A) not(isempty(A)),mtab.type); % find filled rows
mtab  = mtab(idx,:); % keep only filled rows 
mtab  = flip(mtab);

% distinguish in-line and clamp-on (CO)
idx   = find(string(mtab.('bool_clampon'))=='False'); % in-line elements
idxCO = flip(find(string(mtab.('bool_clampon'))=='True')); % clamp-on elemetns

% total number of elements
N       = length(idx);
NCO     = length(idxCO);

% element names
mtab_name = char(mtab.("name"));
moorele   = mtab_name(idx,:);
mooreleCO = mtab_name(idxCO,:);

% modify moorele to correct length
if size(moorele,2)>16
    moorele = moorele(:,1:16);
elseif size(moorele,2)<16
    moorele = char(arrayfun(@(x) x+join(repmat(" ",16-size(moorele,2),1),""),string(moorele)));
end

if size(mooreleCO,2)>16
    mooreleCO = mooreleCO(:,1:16);
elseif size(mooreleCO,2)<16
    mooreleCO = char(arrayfun(@(x) x+join(repmat(" ",16-size(mooreleCO,2),1),""),string(mooreleCO)));
end

% allocate 
[B  ,Cd  ,L  ,Wi  ,D  ,typ  ,ME  ] = deal(zeros(1,N));
[BCO,CdCO,LCO,WiCO,DCO,typCO,ZCO,Jobj,Pobj] = deal(zeros(1,NCO));

% get mooring properties of in-line elements
melem = string(mtab.("name")(idx));
mtype = string(mtab.("type")(idx));
for i=1:N
   if string(mtab.bool_line{idx(i)}) == "True"
       B(i) = mtab.("buoyancy")(idx(i))/mtab.("length")(idx(i));
   else
       B(i) = mtab.("buoyancy")(idx(i));   
   end
   Cd(i) = mtab.("drag")(idx(i));
   Wi(i) = mtab.("width")(idx(i));
   D(i)  = mtab.("diameter")(idx(i));

   mtypei = mtype(i);
   if mtypei == "wires"
        typ(i) = 1;
   elseif mtypei == "chains"
        typ(i) = 2;
   else
        typ(i) = 0;
   end

   mat = mtab.("material")(idx(i));
   if mtypei == "wires"
       if mat==1 % steel
            ME(i)=1.38e11;
       elseif mat==2 % Nylon
            ME(i)=3.45e8;
       elseif mat==3 % Dacron
            ME(i)=8.0e8;
       elseif mat==4 % Polyprop
            ME(i)=3.45e8;
       elseif mat==5 % Polyethy
            ME(i)=6.9e8;
       elseif mat==6 % Kevlar
            ME(i)=6.9e10;
       elseif mat==7 % Aluminum
            ME(i)=7.6e10;
       elseif mat==8 % Dyneema
            ME(i)=1.0e11;
       end
   else
       ME(i) = Inf;
   end

   L(i) = mtab.("length")(idx(i));

end
H = vertcat(L,Wi,D,typ);

% get mooring properties of clamp-on elements
mnameCO = string(mtab.("name")(idxCO));
mtypeCO = string(mtab.("type")(idxCO));
mserCO = string(mtab.("serial")(idxCO));
for i=1:NCO
   BCO(i)   = mtab.("buoyancy")(idxCO(i));
   CdCO(i)  = mtab.("drag")(idxCO(i));
   WiCO(i)  = mtab.("width")(idxCO(i));
   DCO(i)   = mtab.("diameter")(idxCO(i));

   LCO(i)  = mtab.("length")(idxCO(i));
   Jobj(i) = (N+1)-mtab.("num_inline")(idxCO(i));

   section_start = mtab.("inline_bottom")(idxCO(i));
   section_end   = mtab.("inline_top")(idxCO(i));
   section_l     = section_end-section_start;

%    section_h     = str2double(mtab.("SECTION HEIGHT"){idxCO(i)}); % Old language
   section_h     = mtab.("height_along_inline")(idxCO(i));

   Pobj(i) = section_h/section_l;

%    ZCO(i)  = str2double(mtab.("CLAMP MIDDLE HEIGHT")(idxCO(i))); % Old language
   ZCO(i)  = mtab.("height")(idxCO(i));

   mtypeCOi = mtypeCO(i);
   if mtypeCOi == "wires"
        typCO(i) = 1;
   elseif mtypeCOi == "chains"
        typCO(i) = 2;
   else
        typCO(i) = 0;
   end

end
HCO = vertcat(LCO,WiCO,DCO,typCO);

% allocate remainig variables
Iobj = [];
rho  = [];
time = [];
U = [];
V = [];
W = [];
z = [];

%clearvars mtab_name dir file_name i idx idxCO mname mnameCO N NCO psheets mtypeCOi mtypei ptab mat ptabi mtab j L LCO Wi WiCO D DCO typ typCO mtype mtypeCO section_end section_h section_l section_start

%if not(isempty(dir_save))
%    save(dir_save+mname+".mat",'B','BCO','Cd','CdCO','H','HCO','Iobj','Jobj','ME','moorele','mooreleCO','Pobj','rho','time','U','V','W','z','ZCO')
%end

M = struct();
M.B         = B;         
M.BCO       = BCO;       
M.Cd        = Cd;        
M.CdCO      = CdCO;      
M.H         = H;         
M.HCO       = HCO;       
M.Iobj      = Iobj;      
M.Jobj      = Jobj;      
M.ME        = ME;        
M.moorele   = moorele;   
M.mooreleCO = mooreleCO; 
M.Pobj      = Pobj;      
M.rho       = rho;       
M.time      = time;      
M.U         = U;         
M.V         = V;         
M.W         = W;         
M.z         = z;         
M.ZCO       = ZCO;     
M.bonus_mtab=mtab;
M.bonus_ptab=NaN;
end
