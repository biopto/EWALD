function print_plot(handle, filename, dpi, pos)

    set(handle, 'PaperUnits', 'inches', 'PaperPosition', pos/dpi);
    print(handle, '-dpng', ['-r' num2str(dpi)], filename);
    
end
