function alty_video

       hold on 
       
       title(['Aviso: ',num2str(time)])
       imagesc(lon,lat,depth')
       
       saveas( gcf, 'Aviso_',num2str(i), 'jpg' );
       



end