"""
TKPath Compensation Engine

Executes blockchain-verified payments to traditional knowledge holders.
Integrates with OmniPath token infrastructure for transparent distribution.

Trade Secret: Distribution algorithms, payment scheduling, community verification.

Key Innovation: Surcharges from companies opting OUT of TK features fund
compensation to indigenous communities - extractive corporations subsidize justice.
"""

from enum import Enum
from typing import Dict, List, Optional
from pydantic import BaseModel, Field
from datetime import datetime
import uuid
import hashlib

from .attribution import TKContribution, CommunityAttribution


class CompensationStatus(str, Enum):
    """Status of compensation transaction."""
    PENDING = "pending"
    APPROVED = "approved"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"
    DISPUTED = "disputed"


class PaymentMethod(str, Enum):
    """Available payment methods for TK compensation."""
    OMNIPATH_TOKEN = "omnipath_token"      # Native OmniPath blockchain tokens
    STABLECOIN_USDC = "stablecoin_usdc"    # USDC stablecoin
    BANK_TRANSFER = "bank_transfer"         # Traditional bank wire
    COMMUNITY_FUND = "community_fund"       # Pooled community fund
    ESCROW = "escrow"                       # Held in escrow pending verification


class CompensationRecipient(BaseModel):
    """Recipient of compensation payment."""
    recipient_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    community_id: str
    community_name: str
    payment_address: Optional[str] = None  # Blockchain address or account
    payment_method: PaymentMethod = PaymentMethod.COMMUNITY_FUND
    verification_status: str = "verified"  # "pending", "verified", "suspended"
    kyc_completed: bool = False  # Know Your Customer for large payments


class CompensationTransaction(BaseModel):
    """Record of compensation payment transaction."""
    transaction_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    attribution_id: str  # Links to TKContribution
    formulation_id: str
    
    # Payment details
    total_amount: float
    currency: str = "USD"
    recipients: List[CompensationRecipient] = Field(default_factory=list)
    recipient_amounts: Dict[str, float] = Field(default_factory=dict)
    
    # Status tracking
    status: CompensationStatus = CompensationStatus.PENDING
    created_at: datetime = Field(default_factory=datetime.utcnow)
    approved_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    
    # Audit trail
    initiated_by: str = "system"
    approved_by: Optional[str] = None
    blockchain_tx_hash: Optional[str] = None
    verification_hash: Optional[str] = None
    
    # Source tracking
    revenue_source: str = "subscription"  # "subscription", "surcharge", "grant"
    surcharge_funded_percentage: float = 0.0  # % funded by opt-out surcharges
    
    def generate_verification_hash(self) -> str:
        """Generate tamper-evident hash of transaction."""
        data = (
            f"{self.transaction_id}:{self.attribution_id}:{self.total_amount}:"
            f"{self.created_at.isoformat()}"
        )
        for recipient_id, amount in sorted(self.recipient_amounts.items()):
            data += f":{recipient_id}:{amount}"
        self.verification_hash = hashlib.sha256(data.encode()).hexdigest()
        return self.verification_hash


class SurchargePool(BaseModel):
    """Pool of funds collected from TK opt-out surcharges."""
    pool_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    balance: float = 0.0
    total_collected: float = 0.0
    total_distributed: float = 0.0
    last_distribution: Optional[datetime] = None
    
    # Tracking by source
    surcharges_by_company: Dict[str, float] = Field(default_factory=dict)
    distributions_by_community: Dict[str, float] = Field(default_factory=dict)


class TKCompensationEngine:
    """
    Execute blockchain-verified payments to TK holders.
    
    Trade Secret: Distribution algorithms, payment scheduling,
    community verification protocols, surcharge pool management.
    
    Economic Innovation: Companies opting out of TK attribution pay
    surcharges that fund this compensation system - creating a wealth
    transfer from extractive corporations to knowledge holders.
    """
    
    # Distribution parameters (Trade Secret)
    _DISTRIBUTION_SCHEDULE = "monthly"  # Payment frequency
    _MINIMUM_PAYMENT_THRESHOLD = 25.0   # Minimum payout amount
    _PLATFORM_FEE_PERCENTAGE = 5.0      # Platform operational fee
    _ESCROW_HOLD_DAYS = 30              # Days to hold in escrow for disputes
    
    def __init__(self):
        self._transactions: Dict[str, CompensationTransaction] = {}
        self._surcharge_pool = SurchargePool()
        self._pending_distributions: List[CompensationTransaction] = []
    
    def distribute_compensation(
        self,
        revenue: float,
        attribution: TKContribution,
        revenue_source: str = "subscription",
        surcharge_funded: bool = False,
    ) -> CompensationTransaction:
        """
        Execute blockchain-verified payments to TK holders.
        
        Uses OmniPath token infrastructure for transparent distribution.
        
        Args:
            revenue: Total revenue to distribute from
            attribution: TKContribution with community allocations
            revenue_source: Source of funds ("subscription", "surcharge", "grant")
            surcharge_funded: Whether funded by TK opt-out surcharges
            
        Returns:
            CompensationTransaction with payment details
        """
        # Calculate total compensation amount
        compensation_percentage = attribution.suggested_compensation_percentage / 100
        total_compensation = max(
            revenue * compensation_percentage,
            attribution.suggested_minimum_payment
        )
        
        # Apply platform fee
        platform_fee = total_compensation * (self._PLATFORM_FEE_PERCENTAGE / 100)
        distributable_amount = total_compensation - platform_fee
        
        # Create transaction record
        transaction = CompensationTransaction(
            attribution_id=attribution.contribution_id,
            formulation_id=attribution.formulation_id,
            total_amount=total_compensation,
            revenue_source=revenue_source,
            surcharge_funded_percentage=100.0 if surcharge_funded else 0.0,
        )
        
        # Calculate per-community distribution
        total_contribution = sum(
            c.total_contribution_percentage for c in attribution.communities
        )
        
        for community in attribution.communities:
            if total_contribution > 0:
                community_share = (
                    community.total_contribution_percentage / total_contribution
                )
                amount = distributable_amount * community_share
            else:
                # Equal split if no percentage data
                amount = distributable_amount / len(attribution.communities)
            
            # Create recipient record
            recipient = CompensationRecipient(
                community_id=community.community_id,
                community_name=community.community_name,
                payment_address=community.compensation_address,
            )
            
            transaction.recipients.append(recipient)
            transaction.recipient_amounts[recipient.recipient_id] = round(amount, 2)
        
        # Generate verification hash
        transaction.generate_verification_hash()
        
        # Store transaction
        self._transactions[transaction.transaction_id] = transaction
        self._pending_distributions.append(transaction)
        
        return transaction
    
    def process_surcharge_contribution(
        self,
        company_id: str,
        surcharge_amount: float,
    ) -> Dict:
        """
        Add surcharge payment to compensation pool.
        
        When companies opt OUT of TK features, their surcharge payments
        fund compensation to indigenous communities.
        
        Trade Secret: Surcharge allocation algorithm.
        """
        # Add to pool
        self._surcharge_pool.balance += surcharge_amount
        self._surcharge_pool.total_collected += surcharge_amount
        
        # Track by company
        if company_id not in self._surcharge_pool.surcharges_by_company:
            self._surcharge_pool.surcharges_by_company[company_id] = 0.0
        self._surcharge_pool.surcharges_by_company[company_id] += surcharge_amount
        
        return {
            "success": True,
            "company_id": company_id,
            "amount_added": surcharge_amount,
            "pool_balance": self._surcharge_pool.balance,
            "message": (
                f"${surcharge_amount:.2f} surcharge added to TK compensation pool. "
                f"This funds payments to indigenous knowledge holders."
            ),
        }
    
    def execute_monthly_distribution(self) -> List[CompensationTransaction]:
        """
        Execute monthly distribution of surcharge pool to communities.
        
        Trade Secret: Distribution weighting algorithm.
        """
        if self._surcharge_pool.balance < self._MINIMUM_PAYMENT_THRESHOLD:
            return []
        
        # Get all communities that have contributed
        community_weights = self._calculate_community_weights()
        
        if not community_weights:
            return []
        
        distributable = self._surcharge_pool.balance * 0.9  # Keep 10% reserve
        transactions = []
        
        for community_id, weight in community_weights.items():
            amount = distributable * weight
            if amount >= self._MINIMUM_PAYMENT_THRESHOLD:
                tx = CompensationTransaction(
                    attribution_id="monthly_pool_distribution",
                    formulation_id="pool",
                    total_amount=amount,
                    revenue_source="surcharge",
                    surcharge_funded_percentage=100.0,
                )
                
                recipient = CompensationRecipient(
                    community_id=community_id,
                    community_name=community_id,  # Would resolve name in production
                )
                tx.recipients.append(recipient)
                tx.recipient_amounts[recipient.recipient_id] = amount
                tx.status = CompensationStatus.APPROVED
                tx.approved_at = datetime.utcnow()
                
                transactions.append(tx)
                self._transactions[tx.transaction_id] = tx
                
                # Track distribution
                if community_id not in self._surcharge_pool.distributions_by_community:
                    self._surcharge_pool.distributions_by_community[community_id] = 0.0
                self._surcharge_pool.distributions_by_community[community_id] += amount
        
        # Update pool
        total_distributed = sum(tx.total_amount for tx in transactions)
        self._surcharge_pool.balance -= total_distributed
        self._surcharge_pool.total_distributed += total_distributed
        self._surcharge_pool.last_distribution = datetime.utcnow()
        
        return transactions
    
    def _calculate_community_weights(self) -> Dict[str, float]:
        """
        Calculate distribution weights for communities.
        Based on historical attribution frequency.
        Trade Secret: Weighting algorithm.
        """
        # In production, would analyze historical attributions
        # For now, return equal weights for known communities
        from .attribution import TKPathAttribution
        
        attrib = TKPathAttribution()
        communities = attrib.list_known_communities()
        
        if not communities:
            return {}
        
        weight = 1.0 / len(communities)
        return {c["id"]: weight for c in communities}
    
    def approve_transaction(
        self,
        transaction_id: str,
        approver_id: str,
    ) -> Dict:
        """Approve a pending compensation transaction."""
        tx = self._transactions.get(transaction_id)
        if not tx:
            return {"success": False, "error": "Transaction not found"}
        
        if tx.status != CompensationStatus.PENDING:
            return {"success": False, "error": f"Transaction already {tx.status.value}"}
        
        tx.status = CompensationStatus.APPROVED
        tx.approved_by = approver_id
        tx.approved_at = datetime.utcnow()
        
        return {
            "success": True,
            "transaction_id": transaction_id,
            "status": tx.status.value,
            "approved_by": approver_id,
        }
    
    def complete_transaction(
        self,
        transaction_id: str,
        blockchain_tx_hash: Optional[str] = None,
    ) -> Dict:
        """Mark transaction as completed (payment sent)."""
        tx = self._transactions.get(transaction_id)
        if not tx:
            return {"success": False, "error": "Transaction not found"}
        
        if tx.status != CompensationStatus.APPROVED:
            return {"success": False, "error": "Transaction must be approved first"}
        
        tx.status = CompensationStatus.COMPLETED
        tx.completed_at = datetime.utcnow()
        tx.blockchain_tx_hash = blockchain_tx_hash
        
        # Remove from pending
        if tx in self._pending_distributions:
            self._pending_distributions.remove(tx)
        
        return {
            "success": True,
            "transaction_id": transaction_id,
            "status": tx.status.value,
            "completed_at": tx.completed_at.isoformat(),
            "blockchain_tx_hash": blockchain_tx_hash,
        }
    
    def get_transaction(self, transaction_id: str) -> Optional[CompensationTransaction]:
        """Retrieve transaction by ID."""
        return self._transactions.get(transaction_id)
    
    def get_pool_status(self) -> Dict:
        """Get current status of surcharge compensation pool."""
        return {
            "pool_id": self._surcharge_pool.pool_id,
            "current_balance": self._surcharge_pool.balance,
            "total_collected": self._surcharge_pool.total_collected,
            "total_distributed": self._surcharge_pool.total_distributed,
            "pending_distributions": len(self._pending_distributions),
            "last_distribution": (
                self._surcharge_pool.last_distribution.isoformat()
                if self._surcharge_pool.last_distribution else None
            ),
            "contributing_companies": len(self._surcharge_pool.surcharges_by_company),
            "receiving_communities": len(self._surcharge_pool.distributions_by_community),
        }
    
    def get_company_contribution_report(self, company_id: str) -> Dict:
        """
        Get report of surcharge contributions from a specific company.
        Shows how their opt-out payments fund TK compensation.
        """
        amount = self._surcharge_pool.surcharges_by_company.get(company_id, 0.0)
        
        return {
            "company_id": company_id,
            "total_surcharge_contributions": amount,
            "pool_percentage": (
                (amount / self._surcharge_pool.total_collected * 100)
                if self._surcharge_pool.total_collected > 0 else 0.0
            ),
            "message": (
                f"Your TK opt-out surcharges (${amount:.2f}) have been added to the "
                f"Traditional Knowledge Compensation Pool, funding payments to "
                f"indigenous communities whose knowledge informs botanical research."
            ),
        }
    
    def get_audit_trail(
        self,
        company_id: Optional[str] = None,
        community_id: Optional[str] = None,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
    ) -> List[CompensationTransaction]:
        """
        Generate audit trail of compensation transactions.
        For regulatory compliance and transparency reporting.
        """
        transactions = list(self._transactions.values())
        
        if start_date:
            transactions = [t for t in transactions if t.created_at >= start_date]
        if end_date:
            transactions = [t for t in transactions if t.created_at <= end_date]
        if community_id:
            transactions = [
                t for t in transactions
                if any(r.community_id == community_id for r in t.recipients)
            ]
        
        return sorted(transactions, key=lambda t: t.created_at, reverse=True)
