#pragma once

#include "cloth/contact.hpp"
#include "cloth/implicit.hpp"
#include "cloth/mass.hpp"
#include "cloth/state.hpp"
#include "cloth/xpbd.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace cloth {

// The core deliberately has no JSON dependency. A loader (CLI, Python or a
// future service) translates its document into this value before construction.
enum class SceneSolver { implicit_euler, xpbd_distance, xpbd_energy };

struct RectangularSceneSpec {
    std::size_t nx = 2;
    std::size_t ny = 1;
    double width = 1.0;
    double height = 0.7;
};

struct SceneConfigInput {
    RectangularSceneSpec rectangle{};
    // When present, this validated mesh replaces the generated rectangle.
    // Mesh positions are the rest embedding; initial_positions can override
    // them for a deformed starting state.
    std::optional<Mesh> custom_mesh;
    // Empty arrays mean the generated rest embedding and zero velocity. When
    // supplied, both arrays must contain exactly (nx+1)*(ny+1) vertices.
    std::vector<Vec3> initial_positions;
    std::vector<Vec3> initial_velocities;
    std::vector<std::size_t> fixed_vertices;
    Vec3 gravity{0.0, -9.81, 0.0};
    double dt = 0.01;
    std::size_t steps = 8;
    std::uint64_t run_seed = 0;
    Material material{};
    double density = 800.0;
    MassMode mass_mode = MassMode::lumped;
    GroundContactOptions contact{};
    SceneSolver solver = SceneSolver::xpbd_distance;
    NewtonOptions newton_options{};
    XpbdOptions xpbd_options{};
    // Provenance hashes are supplied by the input/asset layer. asset_sha256
    // identifies the mesh/source asset, scene_sha256 the canonical scene
    // document, and state_sha256 the canonical initial state serialization.
    // Empty means unknown; a nonempty field must be 64 lowercase hexadecimal.
    std::string asset_sha256;
    std::string scene_sha256;
    std::string state_sha256;
};

// Validate and canonicalize an input without constructing a solver. The
// canonical fixed list is sorted ascending; duplicate/out-of-range entries,
// non-finite values, invalid dimensions, options or hashes throw
// std::invalid_argument (size arithmetic overflow throws std::length_error).
void validate_scene_config(const SceneConfigInput& input);
bool valid_sha256(const std::string& value) noexcept;
std::vector<std::size_t> canonical_fixed_vertices(
    std::vector<std::size_t> vertices, std::size_t vertex_count);

// Immutable, validated scene data. It owns copies of mesh, state, masses and
// solver options, so later edits to the loader's input cannot change a run.
class SceneConfig {
public:
    explicit SceneConfig(SceneConfigInput input);

    const RectangularSceneSpec& rectangle() const noexcept { return rectangle_; }
    bool has_custom_mesh() const noexcept { return custom_mesh_; }
    const Mesh& mesh() const noexcept { return mesh_; }
    const SimulationState& initial_state() const noexcept { return initial_state_; }
    const std::vector<double>& masses() const noexcept { return masses_; }
    const std::vector<std::size_t>& fixed_vertices() const noexcept { return fixed_vertices_; }
    const Vec3& gravity() const noexcept { return gravity_; }
    double dt() const noexcept { return dt_; }
    std::size_t steps() const noexcept { return steps_; }
    std::uint64_t run_seed() const noexcept { return run_seed_; }
    const Material& material() const noexcept { return material_; }
    double density() const noexcept { return density_; }
    MassMode mass_mode() const noexcept { return mass_mode_; }
    const GroundContactOptions& contact() const noexcept { return contact_; }
    SceneSolver solver() const noexcept { return solver_; }
    const NewtonOptions& newton_options() const noexcept { return newton_options_; }
    const XpbdOptions& xpbd_options() const noexcept { return xpbd_options_; }
    const std::string& asset_sha256() const noexcept { return asset_sha256_; }
    const std::string& scene_sha256() const noexcept { return scene_sha256_; }
    const std::string& state_sha256() const noexcept { return state_sha256_; }

private:
    RectangularSceneSpec rectangle_;
    bool custom_mesh_ = false;
    Mesh mesh_;
    SimulationState initial_state_;
    std::vector<double> masses_;
    std::vector<std::size_t> fixed_vertices_;
    Vec3 gravity_{};
    double dt_ = 0.0;
    std::size_t steps_ = 0;
    std::uint64_t run_seed_ = 0;
    Material material_{};
    double density_ = 0.0;
    MassMode mass_mode_ = MassMode::lumped;
    GroundContactOptions contact_{};
    SceneSolver solver_ = SceneSolver::xpbd_distance;
    NewtonOptions newton_options_{};
    XpbdOptions xpbd_options_{};
    std::string asset_sha256_;
    std::string scene_sha256_;
    std::string state_sha256_;
};

SceneConfigInput default_scene_config_input();
SceneConfig make_default_scene_config();

} // namespace cloth
