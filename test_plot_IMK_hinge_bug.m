plot(hinge_deformation_TH.hinge_1)

% Base shear reactions
[ base_node_reactions ] = fn_xml_read([opensees_dir filesep 'nodal_base_reaction_' 'x' '.xml']);
temp_base_shear = sum(base_node_reactions(1:(end-5),[2,4]),2)';
temp_base_moment = sum(base_node_reactions(1:(end-5),[3,5]),2)';
% Plot Element v hinge forces
close
hold on
plot(element_TH.ele_1.M_TH_1,'DisplayName','Element Force Recorder')
% plot(base_node_reactions(:,5),'DisplayName','Base Shear Recorder x Height')
plot(temp_base_moment,'DisplayName','Base Moment Recorder')
plot(-hinge_force_TH.hinge_1(:,3),'--k','DisplayName','Element Hinge Recorder')
grid on
box on
grid minor
ylabel('base moment (lbs-in)')
legend('location','northwest')

close
hold on
plot(element_TH.ele_1.V_TH_1,'DisplayName','Element Force Recorder')
% plot(-element_TH.ele_1.V_TH_2)
% plot(base_node_reactions(:,2),'DisplayName','Base Shear Recorder')
plot(temp_base_shear,'DisplayName','Base Shear Recorder')
plot(-hinge_force_TH.hinge_1(:,1),'--k','DisplayName','Element Hinge Recorder')
grid on
box on
grid minor
ylabel('base shear (lbs)')
legend('location','northwest')

close
hold on
plot(element_TH.ele_1.M_TH_1)
plot(100*element_TH.ele_1.V_TH_2)
close

% Plot acceleration Compare
hold on
plot(node_TH.node_1_TH.accel_x_abs_TH,'displayName','Base Node')
plot(node_TH.node_2_TH.accel_x_abs_TH,'displayName','Top Node')
grid on
box on
grid minor
ylabel('Abs Acceleration (g)')
legend('location','northwest')
test=5;