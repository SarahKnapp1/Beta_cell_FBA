%input: model, FBA solution of optimzation
%The function gives an overlay in the minerva map based on the fluxs of
%the optimization or generated a overlay file which can be uploed to the
%web.
function minerva_visualization(model,solution_FBA)
    load('minerva.mat')
    minerva.minervaURL = 'http://www.vmh.life/minerva/galaxy.xhtml';
    minerva.map = 'ReconMap-2.01';
    minerva.login = 'user.name';% insert your user name 
    minerva.password = 'minerva.password';%i nsert your password
    minerva.googleLicenseConsent = 'true';
    [serverResponse,curl_str]= buildFluxDistLayoutFix(minerva, model,solution_FBA ,'overlay_Reactions');

    %if did not succssed to generat the overlay becouse the command was to long
    %create a txt file which will conteint the overlay
    if ~contains(serverResponse,'creator')
        overlay_mat=help_overlay(curl_str);
         %write to a txt file
         writematrix(overlay_mat,'Overlay_Reaction.txt','Delimiter','tab')
         disp('generated successfully a overlay file')
    end
end

