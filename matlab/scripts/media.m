function model = media(model, medium)
%% media.m defines the medium constraints we will impose on the genome-scale metabolic model. 
% For every 2g of D-glc, set model.lb to -5. 
% Need to develop formulation to account for other medium mand scale other lbs

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
if ismember({'RPMI'} , medium) % RPMI
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;
    %model.lb(find(ismember(model.rxns, {'EX_gthrd_e'}))) = -5;
    %model.lb(find(ismember(model.rxns, {'EX_btn_e'}))) = -5;
    %model.lb(find(ismember(model.rxns, {'EX_aqcobal_e'}))) = -5;

% If not specified, set to RPMI
elseif ismember({'nan'} , medium) % DMEM w/ low glucose
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.8.html
elseif ismember({'DMEM'} , medium) % DMEM w/ low glucose
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*4.5/2;
    
% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.80.html
elseif ismember({'L15'} , medium) % No glc; some gal  
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = 0; % L15 
    model.lb(find(ismember(model.rxns, {'EX_gal(e)'}))) = -0.9; % LOW GAL

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.83.html
elseif ismember({'McCoy 5A'} , medium) % McCoy 5A
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*3/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.76.html
elseif ismember({'Iscove'} , medium) % IMDM
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*4.5/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.126.html
elseif ismember({'Waymouth'}, medium) % Waymouth
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*5/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.227.html
elseif ismember({'DMEM:F12 (1:1)'}, medium) % 1:1 DMEM and F12
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*3.15/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.64.html
elseif ismember({'HAM F-12'}, medium) % F12
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*1.80/2;

% https://www.sigmaaldrich.com/catalog/search?term=%CE%B1+MEM&interface=All_ZH&N=0&lang=en&region=US&focus=product
% It is different from DMEM, but need to encode glutamine and nucleosides
% into the code...
elseif ismember({'alpha-MEM'}, medium) % Alternative form to MEM
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*4.5/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
elseif ismember({'RPMI w Gln'}, medium) % RPMI ++ L-gln
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.61.html
elseif ismember({'HAM F-10'}, medium) % F10
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*1.1/2;

% No one sells it - calc'd ratio
elseif ismember({'DMEM:RPMI (2:1)'}, medium) % 2:1 DMEM:RPMI
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5/2;

% No one sells it - calc'd ratio
elseif ismember({'MCDB105:M199'}, medium) % MCDB 105: Medium 199
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -(5*0.72*1)/2;

% https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.314.html
elseif ismember({'Williams E Medium'}, medium) % Williams E Medium
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% Same as DMEM:F12 
elseif ismember({'ACL-4'}, medium) % ACL-4
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*3.15/2;  

% No one sells it - calc'd ratio
elseif ismember({'RPMI:F12'}, medium) % 1:1 RPMI and F12
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*1.8/2;   

% No one sells it - calc'd ratio
elseif ismember({'DMEM:Iscove'}, medium) % 1:1 DMEM and Iscove 
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5;

% No one sells it - calc'd ratio  
elseif ismember({'RPMI:Iscove'}, medium) % 1:1 RPMI and Iscove
    model.lb(find(ismember(model.rxns, {'EX_glc(e)'}))) = -5*2.25;
end

end