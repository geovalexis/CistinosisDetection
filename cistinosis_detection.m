
function cistinosis_detection()
    id_pac = input("Cual es el identificador del paciente?: ", 's');
    id_cortes = input("Que cortes quieres analizar?: ");
    pac_path = strcat('D:\OneDrive - URV\Universidad\QuintoAny\Computer Vision\lab\practica\Imatges cistinosis\PAC\', id_pac);
    list_files = dir(pac_path);
    
    for z=1: length(id_cortes)
        images_der = [];
        images_izq = [];
        
        num_corte = num2str(id_cortes(z));
        if id_cortes(z) < 10
            pattern_der = id_pac + "_\d+_\d+_PENTACAM_R_0" + num_corte; 
            pattern_izq = id_pac + "_\d+_\d+_PENTACAM_L_0" + num_corte; 
        else
            pattern_der = id_pac + "_\d+_\d+_PENTACAM_R_" + num_corte;
            pattern_izq = id_pac + "_\d+_\d+_PENTACAM_L_" + num_corte;
        end
        for i=1: length(list_files)
            match = regexp(list_files(i).name, pattern_der);
            if ~isempty(match)
                images_der(end+1) = i;
            end
            match = regexp(list_files(i).name, pattern_izq);
            if ~isempty(match)
                images_izq(end+1) = i;
            end
        end

        results = struct;
        results.izq = struct("fechas", [], "cristales", struct("PS", [], "CS", [], "PP", [], "CP", []), "intensidades", struct("PS", [], "CS", [], "PP", [], "CP", []));
        results.der = struct("fechas", [], "cristales", struct("PS", [], "CS", [], "PP", [], "CP", []), "intensidades", struct("PS", [], "CS", [], "PP", [], "CP", []));
        for i=1: length(images_der)
            image_path_der = string(list_files(images_der(i)).folder) + "\" + string(list_files(images_der(i)).name);
            image_path_izq = string(list_files(images_izq(i)).folder) + "\" + string(list_files(images_izq(i)).name);
            images_path = containers.Map(["izq", "der"], [image_path_izq, image_path_der]);

            for k = keys(images_path)
                % 1. Localización de la cornea
                ima = imread(images_path(k{1}));
                patient_crop = imcrop(ima, [100 150 700 300]);
                patient_grey = rgb2gray(patient_crop);
                patient_ecual = histeq(patient_grey);
                patient_erosion = imerode(patient_ecual, strel('disk', 5));
                patient_erosion_bin = imbinarize(patient_erosion, "adaptive", "sensitivity", 0.01);
                patient_masked = patient_grey .* uint8(patient_erosion_bin);
                patient_bin = imbinarize(patient_masked, 0.1);
                boundingboxes = regionprops(patient_bin, "area", "BoundingBox");   
                [max_area, max_index] = max([boundingboxes.Area]);
                patient_cornea = imcrop(patient_grey, boundingboxes(max_index).BoundingBox);
                % 2. Aplanamiento de la cornea
                patient_cornea_bin = imbinarize(patient_cornea, 0.08);
                for y = 1 : size(patient_cornea, 2)  
                    x_ini = 1;  
                    for x = 1 : size(patient_cornea, 1)
                        if patient_cornea_bin(x, y) == 1 
                            patient_cornea(x_ini, y) = patient_cornea(x,y);
                            if (x_ini  ~= x)  
                                patient_cornea(x, y) = 0; 
                            end
                            x_ini = x_ini + 1;
                        end
                    end
                end
                patient_cornea_bin = imbinarize(patient_cornea, 0.1);
                boundingboxes = regionprops(patient_cornea_bin, "area", "BoundingBox");   
                [max_area, max_index] = max([boundingboxes.Area]);
                patient_cornea_only = imcrop(patient_cornea, boundingboxes(max_index).BoundingBox); 
                % 3. Segmentación en secciones y conteo de cristales y nivel de gris
                [max_y, max_x] = size(patient_cornea_only);
                min_y_superficie = 0;
                min_y_profundo = max_y * 300/550;
                height_superficie = min_y_profundo;
                height_profundo = max_y - min_y_profundo;
                min_x_izq = 0;
                min_x_centro = max_x * 4/12;
                min_x_der = max_x * 8/12;
                width_izq = min_x_centro;
                width_centro = min_x_der - min_x_centro;
                width_der = max_x - min_x_der;

                zones = [min_x_izq, min_y_superficie, width_izq, height_superficie;... %PS
                         min_x_centro, min_y_superficie, width_centro, height_superficie;... %CS   
                         min_x_der, min_y_superficie, width_der, height_superficie;... %PS
                         min_x_izq, min_y_profundo, width_izq, height_profundo;... %PP
                         min_x_centro, min_y_profundo, width_centro, height_profundo;... %CP  
                         min_x_der, min_y_profundo, width_der, height_profundo; %PP
                         ];
                binarized_ima = imbinarize(patient_cornea_only, 0.6);
                for a=1 : length(zones)
                    segmento = imcrop(binarized_ima, zones(a,:)); % Recortamos el segmento
                    [num_pixels, ~] = imhist(segmento); % histograma de cada zona. Num pixels es el recuento de pixeles negros y blancos. 
                    n_cristales(a) = num_pixels(2); % nos quedamos con el nº de blancos
                    cristal_only_segment = imcrop(patient_cornea_only, zones(a,:)) .* uint8(segmento); 
                    [num_pixels_gris, valores_gris] = imhist(cristal_only_segment);
                    suma_valores_grises = 0;
                    for b = 1 : length(valores_gris)
                        suma_valores_grises = suma_valores_grises + valores_gris(b) * num_pixels_gris(b);
                    end
                    if n_cristales(a) ~= 0
                        suma_valores_grises_normalizado(a) = suma_valores_grises / n_cristales(a);
                    else
                        suma_valores_grises_normalizado(a) = 0;
                    end
                end
                results.(k{1}).cristales.PS(end+1) = n_cristales(1) + n_cristales(3);
                results.(k{1}).intensidades.PS(end+1) = suma_valores_grises_normalizado(1) + suma_valores_grises_normalizado(3);
                results.(k{1}).cristales.CS(end+1) = n_cristales(2); 
                results.(k{1}).intensidades.CS(end+1) = suma_valores_grises_normalizado(2);
                results.(k{1}).cristales.PP(end+1) = n_cristales(4) + n_cristales(6);
                results.(k{1}).intensidades.PP(end+1) = suma_valores_grises_normalizado(4) + suma_valores_grises_normalizado(6); 
                results.(k{1}).cristales.CP(end+1) = n_cristales(5); 
                results.(k{1}).intensidades.CP(end+1) = suma_valores_grises_normalizado(5);

                date_str = strsplit(string(list_files(images_der(i)).name), '_');
                results.(k{1}).fechas(end+1) = date_str(2)+date_str(3);
            end           
        end    
        figure(z); plot_results(results); 
    end
end

function plot_results(results)
    subplot(2,2,1);
    dates = datetime(string(results.izq.fechas), 'InputFormat','yyyyMMddhhmmss');
    plot(dates, results.izq.cristales.PS, 'o')
    hold on 
    plot(dates, results.izq.cristales.CS, '*') 
    hold on
    plot(dates, results.izq.cristales.PP, 'x') 
    hold on
    plot(dates, results.izq.cristales.CP, '+') 
    hold off
    legend("PS", "CS", "PP", "CP");
    title("Numeros cristales");
    
    subplot(2,2,3);
    dates = datetime(string(results.izq.fechas), 'InputFormat','yyyyMMddhhmmss');
    plot(dates, results.izq.intensidades.PS, 'o')
    hold on 
    plot(dates, results.izq.intensidades.CS, '*') 
    hold on
    plot(dates, results.izq.intensidades.PP, 'x') 
    hold on
    plot(dates, results.izq.intensidades.CP, '+') 
    hold off
    legend("PS", "CS", "PP", "CP");
    title("Intensidad de grises");
        
    subplot(2,2,2);
    dates = datetime(string(results.der.fechas), 'InputFormat','yyyyMMddhhmmss');
    plot(dates, results.der.cristales.PS, 'o')
    hold on 
    plot(dates, results.der.cristales.CS, '*') 
    hold on
    plot(dates, results.der.cristales.PP, 'x') 
    hold on
    plot(dates, results.der.cristales.CP, '+') 
    hold off
    legend("PS", "CS", "PP", "CP");
    title("Numeros cristales");

    subplot(2,2,4);
    dates = datetime(string(results.der.fechas), 'InputFormat','yyyyMMddhhmmss');
    plot(dates, results.der.intensidades.PS, 'o')
    hold on 
    plot(dates, results.der.intensidades.CS, '*') 
    hold on
    plot(dates, results.der.intensidades.PP, 'x') 
    hold on
    plot(dates, results.der.intensidades.CP, '+') 
    hold off
    legend("PS", "CS", "PP", "CP");
    title("Intensidad de grises");
    
    sgtitle('Ojo Izquierdo                                Ojo Derecho')
end
