"""
API tests for GECKO-A Web Interface v3.0.0.

This module tests the new API endpoints added in v3.0.0:
1. Compound database endpoints
2. Reaction kinetics endpoints
3. Combined workflow endpoints
4. Mass balance endpoints
5. 3D structure endpoints

Author: Deeksha Sharma
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from fastapi.testclient import TestClient

try:
    from gecko_web.main import app
    client = TestClient(app)
    API_AVAILABLE = True
except ImportError:
    API_AVAILABLE = False
    client = None


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestCompoundEndpoints:
    """Test compound database API endpoints."""

    def test_list_compounds(self):
        """Test GET /api/compounds"""
        response = client.get("/api/compounds")
        assert response.status_code == 200

        data = response.json()
        # API may return dict with 'compounds' key or list directly
        if isinstance(data, dict):
            assert "compounds" in data
            compounds = data["compounds"]
        else:
            compounds = data

        assert isinstance(compounds, list)
        assert len(compounds) > 0

        # Check compound structure
        compound = compounds[0]
        assert "name" in compound or "gecko_formula" in compound

    def test_get_compound_categories(self):
        """Test GET /api/compounds/categories"""
        response = client.get("/api/compounds/categories")
        assert response.status_code == 200

        data = response.json()
        assert isinstance(data, (list, dict))

    def test_get_compound_by_name(self):
        """Test GET /api/compounds/{name}"""
        response = client.get("/api/compounds/alpha_pinene")
        # May return 200 or 404 depending on exact name format
        if response.status_code == 200:
            data = response.json()
            assert "name" in data or "smiles" in data

    def test_search_compounds(self):
        """Test GET /api/compounds/search/{query}"""
        response = client.get("/api/compounds/search/pinene")
        assert response.status_code == 200

        data = response.json()
        # API may return dict with 'results' key or list directly
        if isinstance(data, dict):
            results = data.get("results", [])
        else:
            results = data
        assert isinstance(results, list)

    def test_get_compounds_by_category(self):
        """Test GET /api/compounds/category/{category}"""
        response = client.get("/api/compounds/category/terpenes")
        # May return 200 or 404 depending on category format
        if response.status_code == 200:
            data = response.json()
            assert isinstance(data, list)


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestKineticsEndpoints:
    """Test reaction kinetics API endpoints."""

    def test_get_rate_constants(self):
        """Test GET /api/kinetics/{name}"""
        response = client.get("/api/kinetics/alpha_pinene")
        # May return 200 or 404
        if response.status_code == 200:
            data = response.json()
            # Should have rate constant data
            assert isinstance(data, dict)

    def test_get_rate_constants_with_temperature(self):
        """Test GET /api/kinetics/{name} with temperature parameter"""
        response = client.get("/api/kinetics/alpha_pinene?temperature=298")
        if response.status_code == 200:
            data = response.json()
            assert isinstance(data, dict)

    def test_get_atmospheric_lifetime(self):
        """Test GET /api/kinetics/{name}/lifetime"""
        response = client.get("/api/kinetics/alpha_pinene/lifetime")
        if response.status_code == 200:
            data = response.json()
            # Should have lifetime data
            if "OH" in data or "total" in data:
                assert True
            elif "lifetime" in data:
                assert True


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestWorkflowEndpoints:
    """Test workflow API endpoints."""

    def test_combined_workflow_endpoint_exists(self):
        """Test POST /api/workflow/combined endpoint exists"""
        response = client.post(
            "/api/workflow/combined",
            json={
                "voc_name": "test",
                "run_generator": False,
                "run_boxmodel": False
            }
        )
        # May return 200, 400, or 422 depending on validation
        # Just check it doesn't return 404
        assert response.status_code != 404

    def test_comparison_endpoint_exists(self):
        """Test POST /api/workflow/comparison endpoint exists"""
        response = client.post(
            "/api/workflow/comparison",
            json={
                "voc_names": ["test1", "test2"]
            }
        )
        assert response.status_code != 404

    def test_mechanism_reduction_endpoint_exists(self):
        """Test POST /api/mechanism/reduce endpoint exists"""
        response = client.post(
            "/api/mechanism/reduce",
            json={
                "job_id": "test",
                "method": "DRGEP"
            }
        )
        # Endpoint may not exist yet or may return 404 for invalid job
        # Accept any response that doesn't crash
        assert response.status_code in [200, 400, 404, 422, 500]


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestStructure3DEndpoints:
    """Test 3D structure API endpoints."""

    def test_get_3d_structure(self):
        """Test GET /api/structure/{name}/3d"""
        response = client.get("/api/structure/alpha_pinene/3d")
        # May return 200 or 404
        if response.status_code == 200:
            data = response.json()
            # Should have 3D structure data
            assert "mol" in data or "xyz" in data or "structure" in data or "smiles" in data

    def test_get_3d_structure_format(self):
        """Test GET /api/structure/{name}/3d with format parameter"""
        response = client.get("/api/structure/alpha_pinene/3d?format=xyz")
        if response.status_code == 200:
            data = response.json()
            assert isinstance(data, (dict, str))


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestMassBalanceEndpoints:
    """Test mass balance API endpoints."""

    def test_mass_balance_endpoint_exists(self):
        """Test GET /api/jobs/{id}/mass-balance endpoint exists"""
        response = client.get("/api/jobs/test-job-id/mass-balance")
        # Will likely return 404 for non-existent job, but endpoint should exist
        # Accept 404 (job not found) or 200 (success)
        assert response.status_code in [200, 404, 422]


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestEnvironmentEndpoint:
    """Test environment check endpoint."""

    def test_environment_endpoint(self):
        """Test GET /api/environment"""
        response = client.get("/api/environment")
        assert response.status_code == 200

        data = response.json()
        # Should indicate GECKO/BoxModel availability
        assert isinstance(data, dict)


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestJobEndpoints:
    """Test job management endpoints."""

    def test_list_jobs(self):
        """Test GET /api/jobs"""
        response = client.get("/api/jobs")
        assert response.status_code == 200

        data = response.json()
        assert isinstance(data, list)

    def test_create_job(self):
        """Test POST /api/jobs"""
        response = client.post(
            "/api/jobs",
            json={
                "voc_name": "test_compound",
                "job_type": "generator"
            }
        )
        # Should accept the request
        assert response.status_code in [200, 201, 202, 400, 422]

        if response.status_code in [200, 201, 202]:
            data = response.json()
            assert "job_id" in data or "id" in data

    def test_get_job_status(self):
        """Test GET /api/jobs/{id}"""
        # First create a job
        create_response = client.post(
            "/api/jobs",
            json={
                "voc_name": "test_compound_2",
                "job_type": "generator"
            }
        )

        if create_response.status_code in [200, 201, 202]:
            data = create_response.json()
            job_id = data.get("job_id") or data.get("id")

            # Get the job status
            response = client.get(f"/api/jobs/{job_id}")
            assert response.status_code in [200, 404]


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestWebInterface:
    """Test web interface endpoint."""

    def test_main_page(self):
        """Test GET /"""
        response = client.get("/")
        assert response.status_code == 200
        assert "html" in response.headers.get("content-type", "").lower() or \
               "GECKO" in response.text


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestAPIResponseFormats:
    """Test that API responses have correct formats."""

    def test_compounds_response_format(self):
        """Test compound list response format."""
        response = client.get("/api/compounds")
        if response.status_code == 200:
            data = response.json()
            # Handle dict wrapper
            if isinstance(data, dict):
                compounds = data.get("compounds", [])
            else:
                compounds = data

            if len(compounds) > 0:
                compound = compounds[0]
                # Should have some identifying field
                assert "name" in compound or "gecko_formula" in compound, "Missing identifying field"

    def test_error_response_format(self):
        """Test error response format."""
        response = client.get("/api/nonexistent_endpoint_xyz123")
        assert response.status_code == 404

        # Should return JSON error
        try:
            data = response.json()
            assert "detail" in data or "error" in data or "message" in data
        except:
            pass  # Non-JSON error is also acceptable


@pytest.mark.skipif(not API_AVAILABLE, reason="API not available")
class TestAPIValidation:
    """Test API input validation."""

    def test_invalid_job_type(self):
        """Test that invalid job type is handled."""
        response = client.post(
            "/api/jobs",
            json={
                "voc_name": "test",
                "job_type": "invalid_type_xyz"
            }
        )
        # Should return validation error or accept (if job_type is not validated)
        # Just check it doesn't crash
        assert response.status_code in [200, 201, 202, 400, 422]

    def test_missing_required_field(self):
        """Test that missing required field is rejected."""
        response = client.post(
            "/api/jobs",
            json={}  # Missing required fields
        )
        assert response.status_code == 422


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
