
% Finding Acceptable Regions and 5D Combinations

tic 

%% Finding N, ica_fac and va_ratio
% 
% B = [];
% 
% N = 28;
% for ica_fac = 0.6:0.05:0.8;
%     for va_ratio = 1.7:0.05:2;%0.6:0.05:1;
%         for r = 1: size(Re,1)
%             lambda = Re(r, 1);
%             gamma = Re(r,2);
%             [pA_r,Pout_C_r,pV_r]=GenM_ARGL_021018(N, gamma,lambda, ica_fac, va_ratio);
%             for i=1:N
%                 if (pA_r(i,:) >= 32.7*133 && pA_r(i,:) <= 42.7*133) && (pV_r(i,:) >= 8.7*133 && pV_r(i,:) <= 18.7*133)
%                     B = [B; [N lambda gamma va_ratio ica_fac pA_r(i,:)/133 pV_r(i,:)/133]];
%                     Bu = unique(B,'rows');
%                 end
%             end
%         end
%     end
% end


% for N = 27:1:29;
%     for ica_fac = 0.6:0.05:0.8;
%         for va_ratio = 1.7:0.01:2;%0.6:0.05:1;
%             for r = 1: size(Re,1)
%                 lambda = Re(r, 1);
%                 gamma = Re(r,2);
%                 [pA_r,Pout_C_r,pV_r]=GenM_ARGL_021018(N, gamma,lambda, ica_fac, va_ratio);
%                 for i=1:N
%                     if (pA_r(i,:) >= 32.7*133 && pA_r(i,:) <= 42.7*133) && (pV_r(i,:) >= 8.7*133 && pV_r(i,:) <= 18.7*133)
%                         B = [B; [N lambda gamma va_ratio ica_fac pA_r(i,:)/133 pV_r(i,:)/133]];
%                         Bu = unique(B,'rows');
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% % Plots
% for j=1:length(Bu)
%     if (Bu(j,4) == 0.6) && (Bu(j,3) == 0.6)
% %         figure, hold all
%         plot(Bu(j,1),Bu(j,2),'*r','Linewidth',2)
%         hold on, grid on, xlabel('Lambda \lambda'), ylabel('Gamma \gamma'), legend  
%         %         disp(j);
%     end
%     if (Bu(j,4) == 0.7) && (Bu(j,3) == 0.6)
%         plot(Bu(j,1),Bu(j,2),'*b','Linewidth',2)
%         hold on
%     end
%     if (Bu(j,4) == 0.8) && (Bu(j,3) == 0.6)
%         plot(Bu(j,1),Bu(j,2),'*g','Linewidth',2)
%         hold on
%     end
%     if (Bu(j,4) == 0.6) && (Bu(j,3) == 0.7)
%         plot(Bu(j,1),Bu(j,2),'og','Linewidth',2)
%         hold on
%     end
%     if (Bu(j,4) == 0.6) && (Bu(j,3) == 0.8)
%         plot(Bu(j,1),Bu(j,2),'oy','Linewidth',2)
%         hold on
%         legend('avr=0.6, ica-fac=0.6','avr=0.6, ica-fac=0.7','avr=0.6, ica-fac=0.8','avr=0.7, ica-fac=0.6','avr=0.8, ica-fac=0.6')
%     end
% end
% hold off


%% Finding AR data

B = [];
N = 28;
ica_fac = 0.75;
va_ratio= 1.75;

for r = 1: size(Re805,1)
    lambda = Re805(r, 1);
    gamma = Re805(r,2);
    [pA_r,Pout_C_r,pV_r]=GenM_ARGL_101018(N, gamma,lambda, ica_fac, va_ratio);
    for i=1:N
        if (pA_r(i,:) >= 28.5*133 && pA_r(i,:) <= 42.5*133) && (pV_r(i,:) >= -8.6*133 && pV_r(i,:) <= 8.6*133)
%  if (pA_r(i,:) >= 28.5*133 && pA_r(i,:) <= 42.5*133) && (pV_r(i,:)-pV_r(N,:) >= 0*133 && pV_r(i,:)-pV_r(N,:) <= 20*133)
%        if (pA_r(i,:) >= 28.5*133 && pA_r(i,:) <= 42.5*133) && (pV_r(i,:) >= 10.2*133 && pV_r(i,:) <= 24.2*133)
            B = [B; [N lambda gamma va_ratio ica_fac pA_r(i,:)/133 (pV_r(i,:)-pV_r(N,:))/133]];
            Bu = unique(B,'rows');
        end
    end
end



% 
% % Plots
% for j=1:length(Bu)
%     if (Bu(j,4) == 0.6) && (Bu(j,3) == 0.6)
% %         figure, hold all
%         plot(Bu(j,1),Bu(j,2),'*r','Linewidth',2)
%         hold on, grid on, xlabel('Lambda \lambda'), ylabel('Gamma \gamma'), legend  
%         %         disp(j);
%     end
%     if (Bu(j,4) == 0.7) && (Bu(j,3) == 0.6)
%         plot(Bu(j,1),Bu(j,2),'*b','Linewidth',2)
%         hold on
%     end
%     if (Bu(j,4) == 0.8) && (Bu(j,3) == 0.6)
%         plot(Bu(j,1),Bu(j,2),'*g','Linewidth',2)
%         hold on
%     end
%     if (Bu(j,4) == 0.6) && (Bu(j,3) == 0.7)
%         plot(Bu(j,1),Bu(j,2),'og','Linewidth',2)
%         hold on
%     end
%     if (Bu(j,4) == 0.6) && (Bu(j,3) == 0.8)
%         plot(Bu(j,1),Bu(j,2),'oy','Linewidth',2)
%         hold on
%         legend('avr=0.6, ica-fac=0.6','avr=0.6, ica-fac=0.7','avr=0.6, ica-fac=0.8','avr=0.7, ica-fac=0.6','avr=0.8, ica-fac=0.6')
%     end
% end
% hold off

toc

