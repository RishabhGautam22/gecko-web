from fastapi.testclient import TestClient
from gecko_web.main import app

client = TestClient(app)

def test_read_main():
    response = client.get("/")
    assert response.status_code == 200
    assert "GECKO-A Interface" in response.text

def test_create_job():
    response = client.post(
        "/api/jobs",
        json={"voc_name": "test_voc", "job_type": "generator"},
    )
    assert response.status_code == 200
    data = response.json()
    assert "job_id" in data
    assert data["status"] == "queued"

def test_list_jobs():
    # Create a job first
    client.post(
        "/api/jobs",
        json={"voc_name": "test_voc_2", "job_type": "boxmodel"},
    )
    response = client.get("/api/jobs")
    assert response.status_code == 200
    jobs = response.json()
    assert len(jobs) > 0
