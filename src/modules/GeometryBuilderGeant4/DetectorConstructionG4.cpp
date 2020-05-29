/**
 * @file
 * @brief Implements the Geant4 geometry construction process
 * @remarks Code is based on code from Mathieu Benoit
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "DetectorConstructionG4.hpp"

#include <memory>
#include <string>
#include <utility>

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4MultiUnion.hh>
#include <G4NistManager.hh>
#include <G4PVDivision.hh>
#include <G4PVPlacement.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4Sphere.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4VSolid.hh>
#include <G4VisAttributes.hh>
#include "G4Material.hh"

#include "core/geometry/HybridPixelDetectorModel.hpp"
#include "core/module/exceptions.h"
#include "core/utils/log.h"
#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "GeometryConstructionG4.hpp"
#include "Parameterization2DG4.hpp"

using namespace allpix;

DetectorConstructionG4::DetectorConstructionG4(GeometryManager* geo_manager) : geo_manager_(geo_manager) {}

/**
 * @brief Version of std::make_shared that does not delete the pointer
 *
 * This version is needed because some pointers are deleted by Geant4 internally, but they are stored as std::shared_ptr in
 * the framework.
 */
template <typename T, typename... Args> static std::shared_ptr<T> make_shared_no_delete(Args... args) {
    return std::shared_ptr<T>(new T(args...), [](T*) {});
}

void DetectorConstructionG4::build(std::map<std::string, G4Material*> materials_,
                                   const std::shared_ptr<G4LogicalVolume>& world_log) {

    /*
    Build the individual detectors
    */
    std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
    LOG(TRACE) << "Building " << detectors.size() << " device(s)";

    for(auto& detector : detectors) {
        // Get pointer to the model of the detector
        auto model = detector->getModel();

        std::string name = detector->getName();
        LOG(DEBUG) << "Creating Geant4 model for " << name;
        LOG(DEBUG) << " Wrapper dimensions of model: " << Units::display(model->getSize(), {"mm", "um"});
        LOG(TRACE) << " Sensor dimensions: " << model->getSensorSize();
        LOG(TRACE) << " Chip dimensions: " << model->getChipSize();
        LOG(DEBUG) << " Global position and orientation of the detector:";

        // Create the wrapper box and logical volume
        auto wrapper_box = make_shared_no_delete<G4Box>(
            "wrapper_" + name, model->getSize().x() / 2.0, model->getSize().y() / 2.0, model->getSize().z() / 2.0);
        solids_.push_back(wrapper_box);
        auto wrapper_log = make_shared_no_delete<G4LogicalVolume>(
            wrapper_box.get(), materials_["world_material"], "wrapper_" + name + "_log");
        geo_manager_->setExternalObject(name, "wrapper_log", wrapper_log);

        // Get position and orientation
        auto position = detector->getPosition();
        LOG(DEBUG) << " - Position\t\t:\t" << Units::display(position, {"mm", "um"});
        ROOT::Math::Rotation3D orientation = detector->getOrientation();
        std::vector<double> copy_vec(9);
        orientation.GetComponents(copy_vec.begin(), copy_vec.end());
        ROOT::Math::XYZPoint vx, vy, vz;
        orientation.GetComponents(vx, vy, vz);
        auto rotWrapper = std::make_shared<G4RotationMatrix>(copy_vec.data());
        auto wrapperGeoTranslation = toG4Vector(model->getCenter() - model->getGeometricalCenter());
        wrapperGeoTranslation *= *rotWrapper;
        G4ThreeVector posWrapper = toG4Vector(position) - wrapperGeoTranslation;
        geo_manager_->setExternalObject(name, "rotation_matrix", rotWrapper);
        G4Transform3D transform_phys(*rotWrapper, posWrapper);

        G4LogicalVolumeStore* log_volume_store = G4LogicalVolumeStore::GetInstance();
        G4LogicalVolume* world_log_volume = log_volume_store->GetVolume("World_log");

        if(world_log_volume == nullptr) {
            throw ModuleError("Cannot find world volume");
        }

        // Place the wrapper
        auto wrapper_phys = make_shared_no_delete<G4PVPlacement>(
            transform_phys, wrapper_log.get(), "wrapper_" + name + "_phys", world_log.get(), false, 0, true);
        geo_manager_->setExternalObject(name, "wrapper_phys", wrapper_phys);

        LOG(DEBUG) << " Center of the geometry parts relative to the detector wrapper geometric center:";
        /*
                 SENSOR
                 * the sensitive detector is the part that collects the deposits
        */

        // Create the sensor box and logical volume
        auto sensor_box = make_shared_no_delete<G4Box>("sensor_" + name,
                                                       model->getSensorSize().x() / 2.0,
                                                       model->getSensorSize().y() / 2.0,
                                                       model->getSensorSize().z() / 2.0);
        solids_.push_back(sensor_box);

        /* IMPLANTS
         * excise implants from sensor volume and fill them with the implant material
         */
        std::shared_ptr<G4VSolid> sensor_solid = sensor_box;
        if(model->getImplantSize().z() > std::numeric_limits<double>::epsilon()) {
            LOG(TRACE) << "Found implant with non-negligible depth, excising implants from sensor volume.";
            auto implants = model->getImplantSize();

            // Collect all implants in a G4MultiUnion solid to subtract from sensor solid:
            auto implant_union = std::make_shared<G4MultiUnion>();
            solids_.push_back(implant_union);

            for(unsigned int npix_x = 0; npix_x < model->getNPixels().x(); npix_x++) {
                for(unsigned int npix_y = 0; npix_y < model->getNPixels().y(); npix_y++) {
                    // FIXME: We should extend the implant and shift it to avoid fake surfaces
                    auto implant_box = std::make_shared<G4Box>(
                        "implant_box_" + name, implants.x() / 2.0, implants.y() / 2.0, implants.z() / 2.0);
                    solids_.push_back(implant_box);

                    // Calculate transformation for the solid including possible offsets from pixel center
                    auto offset = model->getImplantOffset();
                    G4Transform3D implant_transform(
                        G4RotationMatrix(),
                        G4ThreeVector(
                            -model->getGridSize().x() / 2.0 + (npix_x + 0.5) * model->getPixelSize().x() + offset.x(),
                            -model->getGridSize().y() / 2.0 + (npix_y + 0.5) * model->getPixelSize().y() + offset.y(),
                            (model->getSensorSize().z() - implants.z()) / 2.0));

                    // Add the new solid to the MultiUnion:
                    implant_union->AddNode(*implant_box.get(), implant_transform);
                }
            }

            // Finalize the construction of the multi-union solid:
            implant_union->Voxelize();

            // Obtain implant material from model:
            auto implant_material_iter = materials_.find(model->getImplantMaterial());
            if(implant_material_iter == materials_.end()) {
                throw ModuleError("Cannot construct implants of material '" + model->getImplantMaterial() + "'");
            }
            // Place physical instance of implant extrusion in model (conductor):
            auto implant_log = make_shared_no_delete<G4LogicalVolume>(
                implant_union.get(), implant_material_iter->second, "implants_" + name + "_log");
            geo_manager_->setExternalObject(name, "implants_log", implant_log);

            // Place the implants box
            auto implant_pos = toG4Vector(model->getSensorCenter() - model->getGeometricalCenter());
            LOG(DEBUG) << "  - Implants\t\t:\t" << Units::display(implant_pos, {"mm", "um"});
            auto implant_phys = make_shared_no_delete<G4PVPlacement>(
                nullptr, implant_pos, implant_log.get(), "implants_" + name + "_phys", wrapper_log.get(), false, 0, true);
            geo_manager_->setExternalObject(name, "implants_phys", implant_phys);

            G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(0, 0, 0));
            auto subtraction_solid = std::make_shared<G4SubtractionSolid>(
                "sensor_implant_subtraction_" + name, sensor_box.get(), implant_union.get(), transform);
            solids_.push_back(subtraction_solid);
            sensor_solid = subtraction_solid;
        }

        auto sensor_log =
            make_shared_no_delete<G4LogicalVolume>(sensor_box.get(), materials_["silicon"], "sensor_" + name + "_log");
        geo_manager_->setExternalObject(name, "sensor_log", sensor_log);

        // Place the sensor box
        auto sensor_pos = toG4Vector(model->getSensorCenter() - model->getGeometricalCenter());
        LOG(DEBUG) << "  - Sensor\t\t:\t" << Units::display(sensor_pos, {"mm", "um"});
        auto sensor_phys = make_shared_no_delete<G4PVPlacement>(
            nullptr, sensor_pos, sensor_log.get(), "sensor_" + name + "_phys", wrapper_log.get(), false, 0, true);
        geo_manager_->setExternalObject(name, "sensor_phys", sensor_phys);

        // Create the pixel box and logical volume
        auto pixel_box = make_shared_no_delete<G4Box>("pixel_" + name,
                                                      model->getPixelSize().x() / 2.0,
                                                      model->getPixelSize().y() / 2.0,
                                                      model->getSensorSize().z() / 2.0);
        solids_.push_back(pixel_box);
        auto pixel_log =
            make_shared_no_delete<G4LogicalVolume>(pixel_box.get(), materials_["silicon"], "pixel_" + name + "_log");
        geo_manager_->setExternalObject(name, "pixel_log", pixel_log);

        // Create the parameterization for the pixel grid
        std::shared_ptr<G4VPVParameterisation> pixel_param =
            std::make_shared<Parameterization2DG4>(model->getNPixels().x(),
                                                   model->getPixelSize().x(),
                                                   model->getPixelSize().y(),
                                                   -model->getGridSize().x() / 2.0,
                                                   -model->getGridSize().y() / 2.0,
                                                   0);
        geo_manager_->setExternalObject(name, "pixel_param", pixel_param);

        // WARNING: do not place the actual parameterization, only use it if we need it
        /*
                 CHIP
                 * the chip connected to the bumps bond and the support
        */

        // Construct the chips only if necessary
        if(model->getChipSize().z() > 1e-9) {
            // Create the chip box
            auto chip_box = make_shared_no_delete<G4Box>("chip_" + name,
                                                         model->getChipSize().x() / 2.0,
                                                         model->getChipSize().y() / 2.0,
                                                         model->getChipSize().z() / 2.0);
            solids_.push_back(chip_box);

            // Create the logical volume for the chip
            auto chip_log =
                make_shared_no_delete<G4LogicalVolume>(chip_box.get(), materials_["silicon"], "chip_" + name + "_log");
            geo_manager_->setExternalObject(name, "chip_log", chip_log);

            // Place the chip
            auto chip_pos = toG4Vector(model->getChipCenter() - model->getGeometricalCenter());
            LOG(DEBUG) << "  - Chip\t\t:\t" << Units::display(chip_pos, {"mm", "um"});
            auto chip_phys = make_shared_no_delete<G4PVPlacement>(
                nullptr, chip_pos, chip_log.get(), "chip_" + name + "_phys", wrapper_log.get(), false, 0, true);
            geo_manager_->setExternalObject(name, "chip_phys", chip_phys);
        }

        /*
         * SUPPORT
         * optional layers of support
         */
        auto supports_log = std::make_shared<std::vector<std::shared_ptr<G4LogicalVolume>>>();
        auto supports_phys = std::make_shared<std::vector<std::shared_ptr<G4PVPlacement>>>();
        int support_idx = 0;
        for(auto& layer : model->getSupportLayers()) {
            // Create the box containing the support
            auto support_box = make_shared_no_delete<G4Box>("support_" + name + "_" + std::to_string(support_idx),
                                                            layer.getSize().x() / 2.0,
                                                            layer.getSize().y() / 2.0,
                                                            layer.getSize().z() / 2.0);
            solids_.push_back(support_box);

            std::shared_ptr<G4VSolid> support_solid = support_box;
            if(layer.hasHole()) {
                // NOTE: Double the hole size in the z-direction to ensure no fake surfaces are created
                auto hole_box = make_shared_no_delete<G4Box>("support_" + name + "_hole_" + std::to_string(support_idx),
                                                             layer.getHoleSize().x() / 2.0,
                                                             layer.getHoleSize().y() / 2.0,
                                                             layer.getHoleSize().z());
                solids_.push_back(hole_box);

                G4Transform3D transform(G4RotationMatrix(), toG4Vector(layer.getHoleCenter() - layer.getCenter()));
                auto subtraction_solid = make_shared_no_delete<G4SubtractionSolid>("support_" + name + "_subtraction_" +
                                                                                       std::to_string(support_idx),
                                                                                   support_box.get(),
                                                                                   hole_box.get(),
                                                                                   transform);
                solids_.push_back(subtraction_solid);
                support_solid = subtraction_solid;
            }

            // Create the logical volume for the support
            auto support_material_iter = materials_.find(layer.getMaterial());
            if(support_material_iter == materials_.end()) {
                throw ModuleError("Cannot construct a support layer of material '" + layer.getMaterial() + "'");
            }
            auto support_log =
                make_shared_no_delete<G4LogicalVolume>(support_solid.get(),
                                                       support_material_iter->second,
                                                       "support_" + name + "_log_" + std::to_string(support_idx));
            supports_log->push_back(support_log);

            // Place the support
            auto support_pos = toG4Vector(layer.getCenter() - model->getGeometricalCenter());
            LOG(DEBUG) << "  - Support\t\t:\t" << Units::display(support_pos, {"mm", "um"});
            auto support_phys =
                make_shared_no_delete<G4PVPlacement>(nullptr,
                                                     support_pos,
                                                     support_log.get(),
                                                     "support_" + name + "_phys_" + std::to_string(support_idx),
                                                     wrapper_log.get(),
                                                     false,
                                                     0,
                                                     true);
            supports_phys->push_back(support_phys);

            ++support_idx;
        }
        geo_manager_->setExternalObject(name, "supports_log", supports_log);
        geo_manager_->setExternalObject(name, "supports_phys", supports_phys);

        // Build the bump bonds only for hybrid pixel detectors
        auto hybrid_model = std::dynamic_pointer_cast<HybridPixelDetectorModel>(model);
        if(hybrid_model != nullptr) {

            /*   BUMPS
                    the bump bonds connect the sensor to the readout chip
                */

            // Get parameters from model
            auto bump_height = hybrid_model->getBumpHeight();
            auto bump_sphere_radius = hybrid_model->getBumpSphereRadius();
            auto bump_cylinder_radius = hybrid_model->getBumpCylinderRadius();

            // Create the volume containing the bumps
            auto bump_box = make_shared_no_delete<G4Box>("bump_box_" + name,
                                                         hybrid_model->getSensorSize().x() / 2.0,
                                                         hybrid_model->getSensorSize().y() / 2.0,
                                                         bump_height / 2.);
            solids_.push_back(bump_box);

            // Create the logical wrapper volume
            auto bumps_wrapper_log = make_shared_no_delete<G4LogicalVolume>(
                bump_box.get(), materials_["world_material"], "bumps_wrapper_" + name + "_log");
            geo_manager_->setExternalObject(name, "bumps_wrapper_log", bumps_wrapper_log);

            // Place the general bumps volume
            G4ThreeVector bumps_pos = toG4Vector(hybrid_model->getBumpsCenter() - hybrid_model->getGeometricalCenter());
            LOG(DEBUG) << "  - Bumps\t\t:\t" << Units::display(bumps_pos, {"mm", "um"});
            auto bumps_wrapper_phys = make_shared_no_delete<G4PVPlacement>(nullptr,
                                                                           bumps_pos,
                                                                           bumps_wrapper_log.get(),
                                                                           "bumps_wrapper_" + name + "_phys",
                                                                           wrapper_log.get(),
                                                                           false,
                                                                           0,
                                                                           true);
            geo_manager_->setExternalObject(name, "bumps_wrapper_phys", bumps_wrapper_phys);

            // Create the individual bump solid
            auto bump_sphere = make_shared_no_delete<G4Sphere>(
                "bumps_" + name + "_sphere", 0, bump_sphere_radius, 0, 360 * CLHEP::deg, 0, 360 * CLHEP::deg);
            solids_.push_back(bump_sphere);
            auto bump_tube = make_shared_no_delete<G4Tubs>(
                "bumps_" + name + "_tube", 0., bump_cylinder_radius, bump_height / 2., 0., 360 * CLHEP::deg);
            solids_.push_back(bump_tube);
            auto bump = make_shared_no_delete<G4UnionSolid>("bumps_" + name, bump_sphere.get(), bump_tube.get());
            solids_.push_back(bump);

            // Create the logical volume for the individual bumps
            auto bumps_cell_log =
                make_shared_no_delete<G4LogicalVolume>(bump.get(), materials_["solder"], "bumps_" + name + "_log");
            geo_manager_->setExternalObject(name, "bumps_cell_log", bumps_cell_log);

            // Place the bump bonds grid
            std::shared_ptr<G4VPVParameterisation> bumps_param = std::make_shared<Parameterization2DG4>(
                hybrid_model->getNPixels().x(),
                hybrid_model->getPixelSize().x(),
                hybrid_model->getPixelSize().y(),
                -(hybrid_model->getNPixels().x() * hybrid_model->getPixelSize().x()) / 2.0 +
                    (hybrid_model->getBumpsCenter().x() - hybrid_model->getCenter().x()),
                -(hybrid_model->getNPixels().y() * hybrid_model->getPixelSize().y()) / 2.0 +
                    (hybrid_model->getBumpsCenter().y() - hybrid_model->getCenter().y()),
                0);
            geo_manager_->setExternalObject(name, "bumps_param", bumps_param);

            std::shared_ptr<G4PVParameterised> bumps_param_phys =
                std::make_shared<ParameterisedG4>("bumps_" + name + "_phys",
                                                  bumps_cell_log.get(),
                                                  bumps_wrapper_log.get(),
                                                  kUndefined,
                                                  hybrid_model->getNPixels().x() * hybrid_model->getNPixels().y(),
                                                  bumps_param.get(),
                                                  false);
            geo_manager_->setExternalObject(name, "bumps_param_phys", bumps_param_phys);
        }

        // ALERT: NO COVER LAYER YET

        LOG(TRACE) << " Constructed detector " << detector->getName() << " successfully";
    }
}
