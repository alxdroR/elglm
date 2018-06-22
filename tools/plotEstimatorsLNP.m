function plotEstimatorsLNP(thtaM,bmle,bmele,bapp,mlet,melet)

[nx,tau]=size(thtaM);
thta = fliplr(thtaM); thta = thta(:);

nr = 1; nc = 4;
subplot(nr,nc,1); imagesc(fliplr(thtaM),[min(thtaM(:)) max(thtaM(:))]); axis off;colormap('gray')
title('truth')
subplot(nr,nc,2); imagesc(reshape(bmle,[nx tau]),[min(thtaM(:)) max(thtaM(:))]); axis off;
title('MLE')
text(0,-0.05,['Squared Error = ' num2str(round(100*norm(bmle-thta)^2)/100)],'Units','normalized');

subplot(nr,nc,3); imagesc(reshape(bmele,[nx tau]),[min(thtaM(:)) max(thtaM(:))]); axis off;
title('MELE')
text(0,-0.05,[num2str(round(100*norm(bmele-thta)^2)/100)],'Units','normalized');

subplot(nr,nc,4); imagesc(reshape(bapp,[nx tau]),[min(thtaM(:)) max(thtaM(:))]); axis off;
title('MLE approximation')
text(0,-0.05,{[[num2str(round(100*norm(bapp-thta)^2)/100)]],['computed ' num2str(mlet/melet) ' times faster']},'Units','normalized');


end


